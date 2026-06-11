! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements a crude placeholder for a later implementation of full
!! ice shelf dynamics.
module MOM_ice_shelf_dynamics

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_IS_diag_mediator, only : post_data=>post_IS_data
use MOM_IS_diag_mediator, only : register_diag_field=>register_MOM_IS_diag_field, safe_alloc_ptr
!use MOM_IS_diag_mediator, only : MOM_IS_diag_mediator_init, set_IS_diag_mediator_grid
use MOM_IS_diag_mediator, only : diag_ctrl, time_type, enable_averages, disable_averaging
use MOM_domains, only : MOM_domains_init, clone_MOM_domain
use MOM_domains, only : pass_var, pass_vector, TO_ALL, CGRID_NE, BGRID_NE, AGRID, CORNER, CENTER
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_io, only : file_exists, slasher, MOM_read_data
use MOM_io, only : open_ASCII_file, get_filename_appendix
use MOM_io, only : APPEND_FILE, WRITEONLY_FILE
use MOM_restart, only : register_restart_field, MOM_restart_CS
use MOM_time_manager, only : time_type, get_time, set_time, time_type_to_real, operator(>)
use MOM_time_manager,  only : operator(+), operator(-), operator(*), operator(/)
use MOM_time_manager,  only : operator(/=), operator(<=), operator(>=), operator(<)
use MOM_unit_scaling, only : unit_scale_type, unit_scaling_init
!MJH use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary
use MOM_ice_shelf_state, only : ice_shelf_state
use MOM_coms, only : reproducing_sum, max_across_PEs, min_across_PEs
use MOM_checksums, only : hchksum, qchksum
use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary_channel,initialize_ice_flow_from_file
use MOM_ice_shelf_initialize, only : initialize_ice_shelf_boundary_from_file,initialize_ice_C_basal_friction
use MOM_ice_shelf_initialize, only : initialize_ice_AGlen, initialize_bed_node_from_file
use MOM_ice_shelf_initialize, only : initialize_DG_thickness_from_node_file
implicit none ; private

#include <MOM_memory.h>

public register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn, update_ice_shelf, IS_dynamics_post_data
public ice_time_step_CFL, ice_shelf_dyn_end, change_in_draft, write_ice_shelf_energy
public shelf_advance_front, ice_shelf_min_thickness_calve, calve_to_mask, volume_above_floatation
public reset_DG_to_cellmean_at_cell, reset_DG_to_cellmean_bulk, is_DG_thickness_active
public accumulate_DG_source_rate
public masked_var_grounded

! SSA inner solver flags
integer, parameter :: INNER_CG = 1       !< Conjugate gradient (default)
integer, parameter :: INNER_MINRES = 2   !< MINRES
integer, parameter :: INNER_CR = 3       !< Conjugate residual

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> The control structure for the ice shelf dynamics.
type, public :: ice_shelf_dyn_CS ; private
  real, pointer, dimension(:,:) :: u_shelf => NULL() !< the zonal velocity of the ice shelf/sheet
                                       !! on q-points (B grid) [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: v_shelf => NULL() !< the meridional velocity of the ice shelf/sheet
                                       !! on q-points (B grid) [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: taudx_shelf => NULL() !< the zonal driving stress of the ice shelf/sheet
                                       !! on q-points (C grid) [R L2 T-2 ~> Pa]
  real, pointer, dimension(:,:) :: taudy_shelf => NULL() !< the meridional driving stress of the ice shelf/sheet
                                       !! on q-points (C grid) [R L2 T-2 ~> Pa]
  real, pointer, dimension(:,:) :: sx_shelf => NULL() !< the zonal surface slope of the ice shelf/sheet
                                       !! on q-points (B grid) [nondim]
  real, pointer, dimension(:,:) :: sy_shelf => NULL() !< the meridional surface slope of the ice shelf/sheet
                                       !! on q-points (B grid) [nondim]
  real, pointer, dimension(:,:) :: u_face_mask => NULL() !< mask for velocity boundary conditions on the C-grid
                                       !! u-face - this is because the FEM cares about FACES THAT GET INTEGRATED OVER,
                                       !! not vertices. Will represent boundary conditions on computational boundary
                                       !! (or permanent boundary between fast-moving and near-stagnant ice
                                       !! FOR NOW: 1=interior bdry, 0=no-flow boundary, 2=stress bdry condition,
                                       !! 3=inhomogeneous Dirichlet boundary for u and v, 4=flux boundary: at these
                                       !! faces a flux will be specified which will override velocities; a homogeneous
                                       !! velocity condition will be specified (this seems to give the solver less
                                       !! difficulty)  5=inhomogenous Dirichlet boundary for u only. 6=inhomogenous
                                       !! Dirichlet boundary for v only
  real, pointer, dimension(:,:) :: v_face_mask => NULL()  !< A mask for velocity boundary conditions on the C-grid
                                       !! v-face, with valued defined similarly to u_face_mask, but 5 is Dirichlet for v
                                       !! and 6 is Dirichlet for u
  real, pointer, dimension(:,:) :: u_face_mask_bdry => NULL() !< A duplicate copy of u_face_mask?
  real, pointer, dimension(:,:) :: v_face_mask_bdry => NULL() !< A duplicate copy of v_face_mask?
  real, pointer, dimension(:,:) :: u_flux_bdry_val => NULL() !< The ice volume flux per unit face length into the cell
                                       !! through open boundary u-faces (where u_face_mask=4) [Z L T-1 ~> m2 s-1]
  real, pointer, dimension(:,:) :: v_flux_bdry_val => NULL() !< The ice volume flux per unit face length into the cell
                                       !! through open boundary v-faces (where v_face_mask=4) [Z L T-1 ~> m2 s-1]??
   ! needed where u_face_mask is equal to 4, similarly for v_face_mask
  real, pointer, dimension(:,:) :: umask => NULL()      !< u-mask on the actual degrees of freedom (B grid)
                                       !! 1=normal node, 3=inhomogeneous boundary node,
                                       !!  0 - no flow node (will also get ice-free nodes)
  real, pointer, dimension(:,:) :: vmask => NULL()      !< v-mask on the actual degrees of freedom (B grid)
                                       !! 1=normal node, 3=inhomogeneous boundary node,
                                       !!  0 - no flow node (will also get ice-free nodes)
  real, pointer, dimension(:,:) :: calve_mask => NULL() !< a mask to prevent the ice shelf front from
                                          !! advancing past its initial position (but it may retreat)
  real, pointer, dimension(:,:) :: t_shelf => NULL() !< Vertically integrated temperature in the ice shelf/stream,
                                                     !! on corner-points (B grid) [C ~> degC]
  real, pointer, dimension(:,:) :: tmask => NULL()   !< A mask on tracer points that is 1 where there is ice.
  real, pointer, dimension(:,:,:) :: ice_visc => NULL() !< Area and depth-integrated Glen's law ice viscosity
                                                        !!  (Pa m3 s) in [R L4 Z T-1 ~> kg m2 s-1].
                                                        !!  at either 1 (cell-centered) or 4 quadrature points per cell
  real, pointer, dimension(:,:,:) :: newton_visc_factor => NULL() !< Newton tangent stiffness coefficient:
                                                      !!  (1/n_glen - 1)/2 * ice_visc / eps_e2 at each
                                                      !!  viscosity quadrature point [R L4 Z T ~> kg m2 s]
  real, pointer, dimension(:,:,:) :: newton_str_ux => NULL() !< Longitudinal x-strain-rate ux at each viscosity
                                                      !!  quadrature point for Newton iterations [T-1 ~> s-1]
  real, pointer, dimension(:,:,:) :: newton_str_vy => NULL() !< Longitudinal y-strain-rate vy at each viscosity
                                                      !!  quadrature point for Newton iterations [T-1 ~> s-1]
  real, pointer, dimension(:,:,:) :: newton_str_sh => NULL() !< Engineering shear strain-rate uy+vx at each
                                                      !!  viscosity quadrature point for Newton iterations [T-1 ~> s-1]
  real, pointer, dimension(:,:) :: newton_umid => NULL() !< Cell-averaged zonal velocity u at the current outer
                                                  !! iterate, for Newton basal drag correction [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: newton_vmid => NULL() !< Cell-averaged meridional velocity v at the current
                                                  !! outer iterate, for Newton basal drag correction [L T-1 ~> m s-1]
  real, pointer, dimension(:,:) :: newton_drag_coef => NULL() !< Newton basal drag correction coefficient:
                                         !! 2 * d(basal_trac)/d(|u|^2) * area = d(tau_b_i)/d(u_j) - basal_trac*delta_ij
                                         !! expressed as the u_i*u_j tensor coefficient [R Z T ~> kg m-2 s]
  real, pointer, dimension(:,:) :: AGlen_visc => NULL() !< Ice-stiffness parameter in Glen's law ice viscosity,
                                                      !! often in [Pa-3 s-1] if n_Glen is 3.
  real, pointer, dimension(:,:) :: u_bdry_val => NULL() !< The zonal ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: v_bdry_val => NULL() !< The meridional ice velocity at inflowing boundaries
                                       !! [L yr-1 ~> m yr-1]
  real, pointer, dimension(:,:) :: h_bdry_val => NULL() !< The ice thickness at inflowing boundaries [Z ~> m].
  real, pointer, dimension(:,:) :: t_bdry_val => NULL() !< The ice temperature at inflowing boundaries [C ~> degC].

  real, pointer, dimension(:,:) :: bed_elev => NULL()  !< The bed elevation used for ice dynamics [Z ~> m],
                                                       !! relative to mean sea-level.  This is
                                                       !! the same as G%bathyT+Z_ref, when below sea-level.
                                                       !! Sign convention: positive below sea-level, negative above.
  real, pointer, dimension(:,:) :: bed_node => NULL()  !< Bed elevation reconstructed at B-grid nodes [Z ~> m],
                                                       !! such that the bilinear interpolant over each cell
                                                       !! integrates to the cell-averaged bed_elev. Continuous
                                                       !! across cell boundaries.
  real, pointer, dimension(:,:,:,:) :: h_nodal => NULL() !< DG(1) nodal Q1 thickness per cell at the
                                                       !! 4 corners [Z ~> m]. h_nodal(i,j,a,b) is the value at
                                                       !! cell-local corner (a in {1,2} = west/east, b in {1,2} =
                                                       !! south/north). Each cell owns its own 4 corner values;
                                                       !! jumps across faces are allowed (true DG).
  real, pointer, dimension(:,:,:,:) :: h_flot => NULL() !< Continuous (C0) flotation-gate thickness in the same
                                                       !! per-cell 4-corner layout as h_nodal [Z ~> m]. Each corner
                                                       !! value is the arithmetic mean of the h_nodal values of all
                                                       !! ice cells (hmask 1 or 3) sharing that B-grid node, so
                                                       !! cells sharing a node hold identical values. Used only in
                                                       !! flotation tests (grounded/floating gates, ground_frac,
                                                       !! Coulomb effective pressure) when DG_GL_GATE_CONTINUOUS is
                                                       !! true; never used as a mass or force magnitude. Recomputed
                                                       !! from h_nodal by compute_h_flot at the start of each outer
                                                       !! velocity solve.
  real, pointer, dimension(:,:,:,:) :: Minv_xi => NULL() !< Per-cell 2x2 inverse of the 1D Q1 mass-matrix factor in
                                                       !! the xi direction [L-1 ~> m-1].
  real, pointer, dimension(:,:,:,:) :: Minv_eta => NULL() !< Per-cell 2x2 inverse of the 1D Q1 mass-matrix factor in
                                                       !! the eta direction [L-1 ~> m-1].
  real, pointer, dimension(:,:,:,:) :: cell_mean_w => NULL() !< Per-corner integration weight w(a,b) = int N(a,b)*J
                                                       !! over the reference cell [L2 ~> m2].
  real, pointer, dimension(:,:) :: h_source_rate => NULL() !< Accumulated cell-mean thickness source rate
                                                       !! from external modules (basal melt + surface SMB)
                                                       !! since the last DG advect step [Z T-1 ~> m s-1].
                                                       !! Consumed inside ice_shelf_advect_DG1_nodal: projected
                                                       !! to a continuous Q1 nodal source field and added to
                                                       !! each SSP-RK2 stage RHS. Reset to zero after consumption.
  real, pointer, dimension(:,:) :: h_source_rate_last => NULL() !< Snapshot of h_source_rate as consumed by the
                                                       !! most recent DG advect step [Z T-1 ~> m s-1], kept for
                                                       !! diagnostic posting after h_source_rate has been zeroed.
  real, pointer, dimension(:,:) :: C_basal_friction => NULL()!< Coefficient in sliding law tau_b = C u^(n_basal_fric),
                               !! units of [R L Z T-2 (s m-1)^(n_basal_fric) ~> Pa (s m-1)^(n_basal_fric)]
  real, pointer, dimension(:,:) :: coef_prefactor => NULL() !< Pre-computed area*C_basal_friction*L_T_to_m_s for
                               !! basal friction quadrature evaluation [R L2 Z T-1 ~> kg s-1].
  real, pointer, dimension(:,:) :: fB_elem => NULL()        !< Pre-computed element-level Coulomb fB parameter
                               !! [(T L-1)^CF_PostPeak]; 0 for Weertman.
                               !! Updated each outer iteration by calc_shelf_basal_prefactors.
  real :: alpha_coulomb = 1.0  !< Coulomb prefactor (CF_PostPeak-1)^(CF_PostPeak-1)/CF_PostPeak^CF_PostPeak [nondim]
  real, pointer, dimension(:,:) :: OD_rt => NULL()         !< A running total for calculating OD_av [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac_rt => NULL() !< A running total for calculating ground_frac.
  real, pointer, dimension(:,:) :: OD_av => NULL()         !< The time average open ocean depth [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac => NULL()   !< Fraction of the time a cell is "exposed", i.e. the column
                               !! thickness is below a threshold and interacting with the rock [nondim].  When this
                               !! is 1, the ice-shelf is grounded
  ! float_cond used to be a persistent CS field; it is now derived inline at use sites
  ! from CS%ground_frac (a GL cell is "0 < ground_frac < 1" under GL_regularize=True).
  real, pointer, dimension(:,:,:,:) :: Phi => NULL() !< The gradients of bilinear basis elements at Gaussian
                                                !! 4 quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, pointer, dimension(:,:,:) :: PhiC => NULL()  !< The gradients of bilinear basis elements at 1 cell-centered
                                                !! quadrature point per cell [L-1 ~> m-1].
  real, pointer, dimension(:,:,:) :: Jac => NULL()   !< Jacobian determinant |J_q| = a_q*d_q of the element
                                                !! mapping at each of the 4 Gaussian quadrature points [L2 ~> m2].
                                                !! Equal to G%areaT only for rectangular elements; differs when
                                                !! opposite cell edges have unequal lengths (non-rectangular quads).
  real, pointer, dimension(:,:,:,:,:,:) :: Phisub => NULL() !< Quadrature structure weights at subgridscale
                                                !!  locations for finite element calculations [nondim]
  integer :: OD_rt_counter = 0 !< A counter of the number of contributions to OD_rt.

  real :: velocity_update_time_step !< The time interval over which to update the ice shelf velocity
                    !! using the nonlinear elliptic equation, or 0 to update every timestep [T ~> s].
                    ! DNGoldberg thinks this should be done no more often than about once a day
                    ! (maybe longer) because it will depend on ocean values  that are averaged over
                    ! this time interval, and solving for the equilibrated flow will begin to lose
                    ! meaning if it is done too frequently.
  real :: elapsed_velocity_time  !< The elapsed time since the ice velocities were last updated [T ~> s].

  real :: g_Earth      !< The gravitational acceleration [L2 Z-1 T-2 ~> m s-2].
  real :: density_ice  !< A typical density of ice [R ~> kg m-3].
  real :: Cp_ice       !< The heat capacity of fresh ice [Q C-1 ~> J kg-1 degC-1].

  logical :: advect_shelf !< If true (default), advect ice shelf and evolve thickness
  logical :: reentrant_x !< If true, the domain is zonally reentrant
  logical :: reentrant_y !< If true, the domain is meridionally reentrant
  logical :: alternate_first_direction_IS !< If true, alternate whether the x- or y-direction
                                          !! updates occur first in directionally split parts of the calculation.
  integer :: first_direction_IS !< An integer that indicates which direction is
                                !! to be updated first in directionally split
                                !! parts of the ice sheet calculation (e.g. advection).
  real    :: first_dir_restart_IS = -1.0 !< A real copy of CS%first_direction_IS for use in restart files
  logical :: calc_flux_inout !< If true, calculate the total flux in/out of the domain. This may be required
                             !! for some configurations to calculate flux within a hole in the domain (e.g. at S. Pole)
  integer :: visc_qps !< The number of quadrature points per cell (1 or 4) on which to calculate ice viscosity.
  character(len=40) :: ice_viscosity_compute !< Specifies whether the ice viscosity is computed internally
                                   !! according to Glen's flow law; is constant (for debugging purposes)
                                   !! or using observed strain rates and read from a file
  logical :: shelf_top_slope_bugs !< If true, use directionally inconsistent estimates of the grid
                            !! spacing when calculating the ice shelf surface slope, and underestimate
                            !! slopes near the edge of the ice shelf by a factor of 2.
  logical :: GL_regularize  !< Specifies whether to regularize the floatation condition
                            !! at the grounding line as in Goldberg Holland Schoof 2009
  integer :: n_sub_regularize
                            !< partition of cell over which to integrate for
                            !! interpolated grounding line the (rectangular) is
                            !! divided into nxn equally-sized rectangles, over which
                            !!  basal contribution is integrated (iterative quadrature)
  logical :: GL_couple      !< whether to let the floatation condition be
                            !! determined by ocean column thickness means update_OD_ffrac
                            !! will be called (note: GL_regularize and GL_couple
                            !! should be exclusive)

  real    :: CFL_factor     !< A factor used to limit subcycled advective timestep in uncoupled runs
                            !! i.e. dt <= CFL_factor * min(dx / u) [nondim]

  real :: min_h_shelf !< The minimum ice thickness used during ice dynamics [Z ~> m].
  real :: min_basal_traction !< The minimum basal traction for grounded ice (Pa m-1 s) [R Z T-1 ~> kg m-2 s-1]
  real :: max_surface_slope !< The maximum allowed ice-sheet surface slope (to ignore, set to zero) [nondim]
  real :: min_ice_visc !< The minimum allowed Glen's law ice viscosity (Pa s), in [R L2 T-1 ~> kg m-1 s-1].

  real :: n_glen            !< Nonlinearity exponent in Glen's Law [nondim]
  real :: eps_glen_min      !< Min. strain rate to avoid infinite Glen's law viscosity, [T-1 ~> s-1].
  real :: n_basal_fric      !< Exponent in sliding law tau_b = C u^(m_slide) [nondim]
  logical :: CoulombFriction !< Use Coulomb friction law (Schoof 2005, Gagliardini et al 2007)
  real :: CF_MinN           !< Minimum Coulomb friction effective pressure [R Z L T-2 ~> Pa]
  real :: CF_PostPeak       !< Coulomb friction post peak exponent [nondim]
  real :: CF_Max            !< Coulomb friction maximum coefficient [nondim]
  real :: density_ocean_avg !< A typical ocean density [R ~> kg m-3].  This does not affect ocean
                            !! circulation or thermodynamics.  It is used to estimate the
                            !! gravitational driving force at the shelf front (until we think of
                            !! a better way to do it, but any difference will be negligible).
  real :: thresh_float_col_depth !< The water column depth over which the shelf if considered to be floating
  logical :: moving_shelf_front  !< Specify whether to advance shelf front (and calve).
  logical :: use_DG_thickness     !< If true, use DG(1) representation for ice thickness with
                                  !! unsplit advection scheme and sub-element driving stress quadrature.
  logical :: use_nodal_bed_file   !< If true, read bed elevation at B-grid nodes from NODAL_BED_FILE
                                  !! into CS%bed_node and derive CS%bed_elev by bilinear averaging.
                                  !! Skips reconstruct_bed_to_nodes and the BED_TOPO_FILE read in
                                  !! initialize_ice_flow_from_file. Requires USE_DG_THICKNESS.
  logical :: nodal_positivity     !< If true, apply Liu-style positivity-preserving limiter
                                  !! to the nodal DG(1) thickness corners.
  logical :: dg_hierarchical_lim  !< If true, apply a per-mode hierarchical
                                  !! Zhang-Shu QP-MPP slope limiter with MLP-u2
                                  !! vertex-based bounds to the nodal DG(1)
                                  !! thickness corners between RK stages.
                                  !! Mass is preserved exactly via orthogonalised
                                  !! mode templates against cell_mean_w.
  logical :: dg_lim_isotropic     !< If true (and dg_hierarchical_lim is also true),
                                  !! use isotropic single-phi Zhang-Shu MPP scaling
                                  !! with Park-Kim MLP-u2 vertex-based envelopes
                                  !! instead of the anisotropic per-mode max-product
                                  !! limiter. Single phi per cell scales the full
                                  !! deviation from Hbar at each corner.
  real :: dg_venkat_K             !< Venkatakrishnan-style gradient-proportional
                                  !! slack coefficient for the surface-slope limiter.
                                  !! Adds K * (|dS/dx| + |dS/dy|) * Delta_ref / 2 to
                                  !! the per-corner envelope tolerance, where the
                                  !! cell-mean gradient is estimated by centred
                                  !! differences of S_cell. K = 1 covers the corner-
                                  !! vs-Sbar offset of a globally linear field
                                  !! exactly; K > 1 allows further curvature
                                  !! tolerance. K = 0 disables the term [nondim].
  real, allocatable :: mu_lim_xi(:,:)     !< Cached orthogonalisation offset for the
                                          !! xi-slope mode template, per cell [nondim].
  real, allocatable :: mu_lim_eta(:,:)    !< Cached orthogonalisation offset for the
                                          !! eta-slope mode template, per cell [nondim].
  real, allocatable :: mu_lim_cross(:,:)  !< Cached orthogonalisation offset for the
                                          !! cross mode template, per cell [nondim].
  real, pointer, dimension(:,:) :: dg_lim_phi_xi    => NULL() !< Per-cell xi-slope mode
                                                              !! scaling factor from the
                                                              !! last hierarchical-limiter
                                                              !! call [nondim].
  real, pointer, dimension(:,:) :: dg_lim_phi_eta   => NULL() !< Per-cell eta-slope mode
                                                              !! scaling factor [nondim].
  real, pointer, dimension(:,:) :: dg_lim_phi_cross => NULL() !< Per-cell cross mode
                                                              !! scaling factor [nondim].
  real, pointer, dimension(:,:) :: dg_lim_mass_drift => NULL() !< Per-cell change in
                                                              !! cell-mean thickness from
                                                              !! the limiter [Z ~> m].
                                                              !! Sanity check; should be
                                                              !! ~machine epsilon.
  real, pointer, dimension(:,:) :: dg_lim_phi => NULL()       !< Per-cell limiter-strength
                                                              !! diagnostic [nondim, 0..1].
                                                              !! For isotropic: the unique phi.
                                                              !! For anisotropic: min over the
                                                              !! three mode phis at this cell
                                                              !! (= the most-limiting factor).
  real, pointer, dimension(:,:) :: dg_lim_pk_factor => NULL() !< Park-Kim MLP-u2 smooth-extrema
                                                              !! indicator [nondim, 0 or 1].
                                                              !! 1 where the Park-Kim sign-
                                                              !! consistency check flags the
                                                              !! cell as a smooth extremum and
                                                              !! grants full slack; 0 where the
                                                              !! check rejects and strict MLP-u2
                                                              !! bounds apply.
  real :: dg_art_visc_c_max       !< Peak dimensionless coefficient on the DG(1) artificial-
                                  !! viscosity face flux at fully-shocky faces [nondim].
                                  !! Face coefficient ramps from 0 (smooth) to c_max (shocky)
                                  !! via the relative-jump smoothness indicator
                                  !! r_face = |Delta h_eq| / H_ref: smooth faces
                                  !! receive ~0 damping (preserving optimal mesh convergence
                                  !! on smooth solutions); shocky faces receive c_max scaled
                                  !! by the advective speed + strain-rate floor. c_max = 0
                                  !! disables the viscosity entirely. Typical 0.5-2.0.
  real :: dg_art_visc_advect_coef !< Dimensionless multiplier on the |u_face| advective
                                  !! contribution to u_eff in the DG(1) artificial viscosity
                                  !! [nondim]. Per face, u_eff = advect_coef * |u_face|
                                  !! + strain_coef * eps_e_face * dx_perp. advect_coef = 1
                                  !! (default) preserves the original formulation;
                                  !! advect_coef = 0 drops the |u| term entirely so damping
                                  !! is purely strain-rate-based (yielding a strictly grid-
                                  !! invariant damping timescale at the cost of leaving fast
                                  !! advective regions undamped if eps_e_face is also small).
  real :: dg_art_visc_strain_coef !< Dimensionless coefficient on the velocity-independent
                                  !! strain-rate-scaled diffusivity floor for the DG(1)
                                  !! artificial viscosity [nondim]. Per-face floor velocity is
                                  !! u_floor = strain_coef * eps_e_face * dx_perp, where
                                  !! eps_e_face is the SSA effective strain rate at the face
                                  !! midpoint (2D second invariant). Added to |u_face| in the
                                  !! flux so jumps are still damped at shear-margin / stagnant-
                                  !! interior faces where |u_face| ~ 0 but strain rate is
                                  !! nonzero. strain_coef = 0 (default) recovers the pure
                                  !! velocity-magnitude scaling.
  real :: dg_art_visc_c_min       !< Baseline (smooth-face) DG(1) artificial-viscosity
                                  !! coefficient [nondim]. Nonzero values damp sub-gate jump
                                  !! drift at the cost of first-order dissipation in smooth
                                  !! regions.
  real :: dg_art_visc_r_lo        !< Smoothness-gate lower threshold on the gate ratio
                                  !! r_face [nondim]; faces below this receive only the
                                  !! baseline coefficient.
  real :: dg_art_visc_r_hi        !< Smoothness-gate saturation threshold on the gate ratio
                                  !! r_face [nondim]; faces at or above this receive the
                                  !! full c_max.
  real :: dg_art_visc_kcell       !< Per-cell stability budget for the DG(1) artificial
                                  !! viscosity [nondim]: the sum over a cell's faces of the
                                  !! jump-mode decay rates times dt is held below this by
                                  !! rescaling the cell's face coefficients. SSP-RK2
                                  !! requires < 2.
  logical :: dg_art_visc_wb_harmonic !< If true, the well-balanced equivalent jump uses the
                                  !! harmonic mean of the per-side dh/ds; if false, the
                                  !! arithmetic mean (legacy). The two coincide on uniform-
                                  !! flotation faces; the harmonic mean makes the jump-mode
                                  !! stability rate amplification exactly 2 at every face.
  logical :: dg_art_visc_gate_surface !< If true, the artificial-viscosity smoothness gate
                                  !! ratio is |[s]|/H_ref (surface-cliff fraction of the
                                  !! local mean thickness); if false, |Delta h_eq|/H_ref
                                  !! (legacy, thickness-relative).
  logical :: dg_art_visc_excess_jump !< If true, the DG(1) artificial viscosity damps only
                                  !! the part of the face surface jump in excess of the jump
                                  !! supported by the two cells' mean surfaces; if false,
                                  !! the whole jump (legacy).
  logical :: dg_art_visc_excess_branch_max !< If true, the mean-supported allowance of the
                                  !! excess-jump scheme is the max |[s]| over the admissible
                                  !! flotation-branch assignments of each side's cell mean
                                  !! (admissible = the branch of the mean and the branch of
                                  !! the face trace, both tested at the face-QP bed). Guards
                                  !! against under-built allowances where the face bed is a
                                  !! local extremum unrepresentative of the cell interiors
                                  !! (e.g. a grounding line on a sill); identical wherever
                                  !! mean and trace agree on the branch.
  real :: dg_slow_idle_u_tiny     !< Stagnant-jump diagnostic threshold on face-speed
                                  !! magnitude [L T-1]. A face is flagged if |u_face| <
                                  !! this value AND eps_e_face < eps_tiny AND |Delta h_eq|
                                  !! > s_tol; persistent flags mark regions where neither
                                  !! the |u| nor the strain-rate term damps the jump mode.
  real :: dg_slow_idle_eps_tiny   !< Stagnant-jump diagnostic threshold on face strain rate
                                  !! [T-1]. See dg_slow_idle_u_tiny.
  real :: dg_slow_idle_s_tol      !< Stagnant-jump diagnostic threshold on |Delta h_eq|
                                  !! (well-balanced thickness-equivalent surface jump) [Z].
                                  !! See dg_slow_idle_u_tiny.
  logical :: dg_driving_stress_IBP !< If true, the DG(1) driving stress uses the
                                  !! integration-by-parts weak form with central P*
                                  !! (= 1/2(P_loc + P_ngh)) at interior faces. If false
                                  !! (default), use a strong-form collocation that evaluates
                                  !! rho*g*h*grad(s) directly at Gauss points using the Q1
                                  !! nodal basis, with an optional scale-aware face flux
                                  !! gated by dg_face_flux_K_thresh.
  real :: dg_face_flux_K_thresh   !< Controls the strong-form driving stress's optional
                                  !! face flux at interior faces (between two hmask=1
                                  !! cells). <0 (default) disables it entirely (pure
                                  !! strong form). =0 forces s=1 (full central-IBP face
                                  !! flux: equivalent at the SSA node assembly to the
                                  !! IBP routine). >0 enables a Venkatakrishnan-style
                                  !! blend s = r^2 / (r^2 + K^2), r = |[h]|/h_avg, that
                                  !! ramps from zero (smooth faces) to the central-IBP
                                  !! flux (large jumps). ~0.05-0.2 engages on real h
                                  !! discontinuities; >>1 is effectively disabled
                                  !! (r << 1 in typical flows) [nondim].
  logical :: dg_gl_gate_continuous !< If true, evaluate all grounded/floating decisions
                                  !! (basal-friction gates, ground_frac, Coulomb effective
                                  !! pressure, and the strong-form driving-stress branch
                                  !! selector) on the continuous node-averaged thickness
                                  !! field h_flot instead of each cell's own DG h_nodal.
                                  !! Prevents the flotation crossing from hiding inside an
                                  !! inter-cell thickness jump, which otherwise lets the
                                  !! grounding line lock at cell faces with ground_frac
                                  !! exactly 0/1 and defeats the subgrid GL scheme. Force
                                  !! magnitudes (rho*g*h, grad h, face-jump terms) always
                                  !! use the true DG h_nodal. Default false.
  logical :: dg_gl_gate_driving_stress !< If true (and dg_gl_gate_continuous is true), the
                                  !! strong-form driving-stress flotation branches (volume,
                                  !! subgrid, and face Dirac term) also use h_flot. If false
                                  !! (default), the driving stress keeps per-side own-h
                                  !! branching, which is equivalent to the continuous
                                  !! s = max(h - b, (1-r)*h) and hence yields forces that are
                                  !! continuous in the state. h_flot gating of the face Dirac
                                  !! term makes [s] switch discontinuously between [h] and
                                  !! (1-r)*[h] at gate sign flips, creating a strong restoring
                                  !! force that pins the steady grounding line at cell faces.
  real :: dg_gl_gate_deficit_scale !< Regularization scale s for inverse-flotation-deficit
                                  !! weighting of the h_flot node average [Z ~> m]. <= 0
                                  !! (default 0) gives the plain arithmetic corner mean; > 0
                                  !! weights each corner by 1/(|r*h - bed| + s) so the side
                                  !! nearer flotation dominates, preventing a large one-sided
                                  !! thickness jump at the GL face from dragging the gate's
                                  !! flotation crossing far into the lighter cell. Small s
                                  !! approaches pinning straddling faces at flotation (dead
                                  !! band; avoid); large s approaches the arithmetic mean.
  logical :: calve_to_mask       !< If true, calve off the ice shelf when it passes the edge of a mask.
  real :: min_thickness_simple_calve !< min. ice shelf thickness criteria for calving [Z ~> m].
  real :: T_shelf_missing   !< An ice shelf temperature to use where there is no ice shelf [C ~> degC]
  real :: cg_tolerance !< For Picard iterations, the tolerance in the CG solver, relative to initial residual, that
                       !! determines when to stop the conjugate gradient iterations [nondim].
  real :: cg_newton_tolerance !< For inexact Newton iterations, the initial tolerance in the CG solver, relative to
                              !!  initial residual, that determines when to stop the CG iterations [nondim].
  real :: cg_tol_current !< Working CG tolerance for the current inner solve [nondim].
  real :: nonlinear_tolerance !< The fractional nonlinear tolerance, relative to the initial error,
                              !! that sets when to stop the iterative velocity solver [nondim]
  real :: newton_after_tolerance !< The fractional nonlinear tolerance, relative to the initial error, at
                                 !! which to switch from Picard to Newton iterations in the velocity solver
                                 !! If set to <= 0, no Picard [nondim]
  logical :: newton_adapt_cg_tol !< Use an adaptive CG tolerance during Newton iterations
  real :: ew_gamma !< Gamma in Eisenstat-Walker adaptive Newton tolerance [nondim].
  real :: ew_alpha !< Alpha in Eisenstat-Walker adaptive Newton tolerance [nondim].
  integer :: ew_safety !< Safeguard Eisenstat-Walker using:
                    !!(0) no safeguard, (1) EW choice 2 threshold or (2) PETSc option 3 (Chacon 2008)
  real :: ew_1_thres !< Threshold for Eisenstat-Walker version 1 [nondim]
  real :: ew_eta_max !< Maximum allowed Eisenstat-Walker eta [nondim]
  integer :: cg_max_iterations !< The maximum number of iterations that can be used in the CG solver
  integer :: nonlin_solve_err_mode  !< 1: exit based on nonlin residual | F | / | F_0 | where | | is infty-norm
                    !! 2: exit based on "fixed point" metric (|u - u_last| / |u| < tol) where | | is infty-norm
                    !! 3: exit based on change of solution norm 2*abs(|u|-|u_last|)/(|u|+|u_last|) where | | is L2-norm
                    !! 4: exit based on nonlin residual  | F | / | F_0 | where | | is L2-norm
                    !! 5: exit based on relative residual | F | / | tau | where | | is L2-norm
  logical :: ssa_add_rel_resid !< Nonlinear error in velocity solve will also depend on the
                                   !! L2 residual norm relative to RHS norm
  real :: rr_nonlinear_tolerance !< If ssa_add_rel_resid, the additional nonlin tolerance in the iterative
                    !! velocity solve used for the relative residual [nondim]
  ! for write_ice_shelf_energy
  type(time_type) :: energysavedays            !< The interval between writing the energies
                                               !! and other integral quantities of the run.
  type(time_type) :: energysavedays_geometric  !< The starting interval for computing a geometric
                                               !! progression of time deltas between calls to
                                               !! write_energy. This interval will increase by a factor of 2.
                                               !! after each call to write_energy.
  logical         :: energysave_geometric      !< Logical to control whether calls to write_energy should
                                               !! follow a geometric progression
  type(time_type) :: write_energy_time         !< The next time to write to the energy file.
  type(time_type) :: geometric_end_time        !< Time at which to stop the geometric progression
                                               !! of calls to write_energy and revert to the standard
                                               !! energysavedays interval
  real    :: timeunit           !< The length of the units for the time axis and certain input parameters
                                !! including ENERGYSAVEDAYS [s].
  type(time_type) :: Start_time !< The start time of the simulation.
                                ! Start_time is set in MOM_initialization.F90
  integer :: prev_IS_energy_calls = 0 !< The number of times write_ice_shelf_energy has been called.
  integer :: IS_fileenergy_ascii   !< The unit number of the ascii version of the energy file.
  character(len=200) :: IS_energyfile  !< The name of the ice sheet energy file with path.

  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1

  logical :: debug                !< If true, write verbose checksums for debugging purposes
                                  !! and use reproducible sums
  logical :: doing_newton = .false. !< If true, the outer iteration is using Newton (tangent) linearization
                                    !! instead of Picard (secant) linearization for the ice viscosity
  integer :: inner_solver !< The inner linear solver: INNER_CG (1),INNER_MINRES (2), or INNER_CR (3)
  logical :: cg_halo_shrink = .true. !< If true, CG uses halo-shrinking to defer pass_vector calls;
                                     !! if false, uses fixed CG_action range with 1 pass_vector per iteration
  logical :: module_is_initialized = .false. !< True if this module has been initialized.

  !>@{ Diagnostic handles
  integer :: id_u_shelf = -1, id_v_shelf = -1, id_shelf_speed, id_t_shelf = -1, &
             id_taudx_shelf = -1, id_taudy_shelf = -1, id_taud_shelf = -1, id_bed_elev = -1, &
             id_ground_frac = -1, id_col_thick = -1, id_OD_av = -1, &
             id_u_mask = -1, id_v_mask = -1, id_ufb_mask =-1, id_vfb_mask = -1, id_t_mask = -1, &
             id_sx_shelf = -1, id_sy_shelf = -1, id_surf_slope_mag_shelf, &
             id_duHdx = -1, id_dvHdy = -1, id_fluxdiv = -1, &
             id_strainrate_xx = -1, id_strainrate_yy = -1, id_strainrate_xy = -1, &
             id_pstrainrate_1 = -1, id_pstrainrate_2, &
             id_devstress_xx = -1, id_devstress_yy = -1, id_devstress_xy = -1, &
             id_pdevstress_1 = -1, id_pdevstress_2 = -1

  !>@}
  ! ids for outputting intermediate thickness in advection subroutine (debugging)
  !>@{ Diagnostic handles for debugging
  integer :: id_h_after_uflux = -1, id_h_after_vflux = -1, id_h_after_adv = -1, &
             id_visc_shelf = -1, id_taub = -1, &
             id_bed_node = -1, &
             id_h_nodal_SW = -1, id_h_nodal_SE = -1, id_h_nodal_NW = -1, id_h_nodal_NE = -1, &
             id_h_jump_node = -1, id_h_jump_node_rel = -1, &
             id_h_node_max = -1, id_h_node_min = -1, &
             id_h_jump_envelope = -1, id_h_jump_envelope_rel = -1, &
             id_h_overshoot_node = -1, id_h_overshoot_node_rel = -1, &
             id_s_overshoot_node = -1, id_s_overshoot_node_rel = -1, id_h_source_rate = -1, &
             id_dg_lim_phi_xi = -1, id_dg_lim_phi_eta = -1, id_dg_lim_phi_cross = -1, &
             id_dg_lim_mass_drift = -1, id_dg_lim_phi = -1, id_dg_lim_pk_factor = -1, &
             id_phi_x_FV = -1, id_phi_y_FV = -1, &
             id_dg_art_visc_coef_u = -1, id_dg_art_visc_coef_v = -1, &
             id_dg_art_visc_nu_u = -1, id_dg_art_visc_nu_v = -1, &
             id_dg_art_visc_excess_frac_u = -1, id_dg_art_visc_excess_frac_v = -1, &
             id_dg_art_visc_allow_u = -1, id_dg_art_visc_allow_v = -1, &
             id_dg_art_visc_cell_scale = -1, &
             id_dg_slow_idle_face_u = -1, id_dg_slow_idle_face_v = -1, &
             id_h_jump_face_u = -1, id_h_jump_face_v = -1, &
             id_s_jump_face_u = -1, id_s_jump_face_v = -1, &
             id_s_jump_face_u_rel = -1, id_s_jump_face_v_rel = -1, &
             id_h_jump_face_u_signed = -1, id_h_jump_face_v_signed = -1, &
             id_s_jump_face_u_signed = -1, id_s_jump_face_v_signed = -1, &
             id_un_face_u = -1, id_un_face_v = -1, &
             id_dg_eps_face_u = -1, id_dg_eps_face_v = -1
  real, pointer, dimension(:,:) :: dg_art_visc_coef_u => NULL() !< Per-face DG(1) artificial-
                                                       !! viscosity coefficient on u-faces
                                                       !! [nondim] (post per-cell CFL
                                                       !! scaling), from the last spatial-
                                                       !! operator call. Diagnostic only.
  real, pointer, dimension(:,:) :: dg_art_visc_coef_v => NULL() !< As dg_art_visc_coef_u but
                                                       !! on v-faces [nondim].
  real, pointer, dimension(:,:) :: dg_art_visc_nu_u => NULL() !< Per-face effective DG(1)
                                                       !! artificial viscosity on u-faces
                                                       !! [m2 s-1], = c_face*u_eff*dx_perp
                                                       !! (post per-cell CFL scaling).
                                                       !! Comparable to a physical diffusivity:
                                                       !! the face flux equals nu * [h_eq] /
                                                       !! dx_perp * dy_face. Distinguishes
                                                       !! gate-active-but-quiet faces from
                                                       !! gate-active-and-damping faces.
  real, pointer, dimension(:,:) :: dg_art_visc_nu_v => NULL() !< As dg_art_visc_nu_u but on
                                                       !! v-faces [m2 s-1].
  real, pointer, dimension(:,:) :: dg_art_visc_excess_frac_u => NULL() !< Fraction of the
                                                       !! well-balanced surface jump treated
                                                       !! as excess by the EXCESS_JUMP
                                                       !! allowance on u-faces [nondim],
                                                       !! max over the 2 face QPs of
                                                       !! |ds_use|/|ds_qp|. Diagnostic only.
  real, pointer, dimension(:,:) :: dg_art_visc_excess_frac_v => NULL() !< As
                                                       !! dg_art_visc_excess_frac_u but on
                                                       !! v-faces [nondim].
  real, pointer, dimension(:,:) :: dg_art_visc_allow_u => NULL() !< Mean-supported surface-
                                                       !! jump allowance |ds_bar| on u-faces
                                                       !! [Z ~> m], max over the 2 face QPs.
                                                       !! Diagnostic only.
  real, pointer, dimension(:,:) :: dg_art_visc_allow_v => NULL() !< As dg_art_visc_allow_u
                                                       !! but on v-faces [Z ~> m].
  real, pointer, dimension(:,:) :: dg_art_visc_cell_scale => NULL() !< Per-cell DG(1)
                                                       !! artificial-viscosity cap throttle
                                                       !! factor [nondim]; 1 where the
                                                       !! per-cell stability budget is slack,
                                                       !! kcell/(S_K*dt) < 1 where the cap
                                                       !! rescales the cell's face
                                                       !! coefficients.
  real, pointer, dimension(:,:) :: dg_slow_idle_face_u => NULL() !< Stagnant-jump indicator on
                                                       !! u-faces [nondim, 0 or 1]. 1 where
                                                       !! |u_face| < u_tiny AND eps_e_face
                                                       !! < eps_tiny AND |Delta h_eq| > s_tol
                                                       !! at the last spatial-operator call,
                                                       !! flagging faces where neither the
                                                       !! advective nor the strain-rate
                                                       !! damping channel acts. Persistent
                                                       !! non-zero values indicate that a
                                                       !! constant velocity-floor B1 may
                                                       !! be warranted.
  real, pointer, dimension(:,:) :: dg_slow_idle_face_v => NULL() !< As dg_slow_idle_face_u but
                                                       !! on v-faces [nondim, 0 or 1].
  real, pointer, dimension(:,:) :: phi_x_FV => NULL() !< Van Leer slope-limiter factor at each u-face
                                                       !! from ice_shelf_advect_thickness_x [nondim],
                                                       !! in [0,2]. Faces where the limiter branch was
                                                       !! not taken (incomplete stencil, inactive face)
                                                       !! report 1.0.
  real, pointer, dimension(:,:) :: phi_y_FV => NULL() !< Van Leer slope-limiter factor at each v-face
                                                       !! from ice_shelf_advect_thickness_y [nondim].
  !>@}
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to control diagnostic output.

end type ice_shelf_dyn_CS

!> A container for loop bounds
type :: loop_bounds_type ; private
  integer :: ish !< Starting i-index of the computational domain [nondim]
  integer :: ieh !< Ending i-index of the computational domain [nondim]
  integer :: jsh !< Starting j-index of the computational domain [nondim]
  integer :: jeh !< Ending j-index of the computational domain [nondim]
end type loop_bounds_type

contains

!> used for flux limiting in advective subroutines Van Leer limiter (source: Wikipedia)
!! The return value is between 0 and 2 [nondim].
function slope_limiter(num, denom)
  real, intent(in)    :: num   !< The numerator of the ratio used in the Van Leer slope limiter
  real, intent(in)    :: denom !< The denominator of the ratio used in the Van Leer slope limiter
  real :: slope_limiter ! The slope limiter value, between 0 and 2 [nondim].
  real :: r  ! The ratio of num/denom [nondim]

  if (denom == 0) then
    slope_limiter = 0
  elseif (num*denom <= 0) then
    slope_limiter = 0
  else
    r = num/denom
    slope_limiter = (r+abs(r))/(1+abs(r))
  endif

end function slope_limiter

!> Calculate area of quadrilateral.
function quad_area (X, Y)
  real, dimension(4), intent(in) :: X !< The x-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(4), intent(in) :: Y !< The y-positions of the vertices of the quadrilateral [L ~> m].
  real :: quad_area ! Computed area [L2 ~> m2]
  real :: p2, q2, a2, c2, b2, d2

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

  p2 = ( ((X(4)-X(1))**2) + ((Y(4)-Y(1))**2) ) ; q2 = ( ((X(3)-X(2))**2) + ((Y(3)-Y(2))**2) )
  a2 = ( ((X(3)-X(4))**2) + ((Y(3)-Y(4))**2) ) ; c2 = ( ((X(1)-X(2))**2) + ((Y(1)-Y(2))**2) )
  b2 = ( ((X(2)-X(4))**2) + ((Y(2)-Y(4))**2) ) ; d2 = ( ((X(3)-X(1))**2) + ((Y(3)-Y(1))**2) )
  quad_area = .25 * sqrt(4*P2*Q2-(B2+D2-A2-C2)**2)

end function quad_area

!> This subroutine is used to register any fields related to the ice shelf
!! dynamics that should be written to or read from the restart file.
subroutine register_ice_shelf_dyn_restarts(G, US, param_file, CS, restart_CS)
  type(ocean_grid_type),  intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(unit_scale_type),  intent(in)    :: US   !< A structure containing unit conversion factors
  type(param_file_type),  intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(ice_shelf_dyn_CS), pointer       :: CS !< A pointer to the ice shelf dynamics control structure
  type(MOM_restart_CS),   intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  real :: T_shelf_missing ! An ice shelf temperature to use where there is no ice shelf [C ~> degC]
  logical :: shelf_mass_is_dynamic, override_shelf_movement, active_shelf_dynamics
  character(len=40)  :: mdl = "MOM_ice_shelf_dyn"  ! This module's name.
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (associated(CS)) then
    call MOM_error(FATAL, "MOM_ice_shelf_dyn.F90, register_ice_shelf_dyn_restarts: "// &
                          "called with an associated control structure.")
    return
  endif
  allocate(CS)

  override_shelf_movement = .false. ; active_shelf_dynamics = .false.
  call get_param(param_file, mdl, "DYNAMIC_SHELF_MASS", shelf_mass_is_dynamic, &
                 "If true, the ice sheet mass can evolve with time.", &
                 default=.false., do_not_log=.true.)
  if (shelf_mass_is_dynamic) then
    call get_param(param_file, mdl, "OVERRIDE_SHELF_MOVEMENT", override_shelf_movement, &
                 "If true, user provided code specifies the ice-shelf "//&
                 "movement instead of the dynamic ice model.", default=.false., do_not_log=.true.)
    active_shelf_dynamics = .not.override_shelf_movement
  endif

  if (active_shelf_dynamics) then
    call get_param(param_file, mdl, "MISSING_SHELF_TEMPERATURE", T_shelf_missing, &
                 "An ice shelf temperature to use where there is no ice shelf.",&
                 units="degC", default=-10.0, scale=US%degC_to_C, do_not_log=.true.)

    call get_param(param_file, mdl, "NUMBER_OF_ICE_VISCOSITY_QUADRATURE_POINTS", CS%visc_qps, &
                 "Number of ice viscosity quadrature points. Either 1 (cell-centered) for 4", &
                  units="none", default=1)
    if (CS%visc_qps/=1 .and. CS%visc_qps/=4) call MOM_error (FATAL, &
      "NUMBER OF ICE_VISCOSITY_QUADRATURE_POINTS must be 1 or 4")

    call get_param(param_file, mdl, "FIRST_DIRECTION_IS", CS%first_direction_IS, &
                 "An integer that indicates which direction goes first "//&
                 "in parts of the code that use directionally split "//&
                 "updates (e.g. advection), with even numbers (or 0) used for x- first "//&
                 "and odd numbers used for y-first.", default=0)
    call get_param(param_file, mdl, "ALTERNATE_FIRST_DIRECTION_IS", CS%alternate_first_direction_IS, &
                 "If true, after every advection call, alternate whether the x- or y- "//&
                 "direction advection updates occur first. "//&
                 "If this is true, FIRST_DIRECTION applies at the start of a new run or if "//&
                 "the next first direction can not be found in the restart file.", default=.false.)
    call get_param(param_file, mdl, "CALC_FLUX_INOUT", CS%calc_flux_inout, &
                 "If true, during every advection call, calculate and output the total flux in/out " //&
                 "of the domain (e.g. at S. Pole)", default=.false.)

    allocate(CS%u_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%v_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%t_shelf(isd:ied,jsd:jed), source=T_shelf_missing) ! [C ~> degC]
    allocate(CS%ice_visc(isd:ied,jsd:jed,CS%visc_qps), source=0.0)
    allocate(CS%newton_visc_factor(isd:ied,jsd:jed,CS%visc_qps), source=0.0)
    allocate(CS%newton_str_ux(isd:ied,jsd:jed,CS%visc_qps), source=0.0)
    allocate(CS%newton_str_vy(isd:ied,jsd:jed,CS%visc_qps), source=0.0)
    allocate(CS%newton_str_sh(isd:ied,jsd:jed,CS%visc_qps), source=0.0)
    allocate(CS%newton_umid(isd:ied,jsd:jed), source=0.0)
    allocate(CS%newton_vmid(isd:ied,jsd:jed), source=0.0)
    allocate(CS%newton_drag_coef(isd:ied,jsd:jed), source=0.0)
    allocate(CS%AGlen_visc(isd:ied,jsd:jed), source=2.261e-25) ! [Pa-3 s-1]
    allocate(CS%C_basal_friction(isd:ied,jsd:jed), source=5.0e10*US%Pa_to_RLZ_T2)
             ! Units of [R L Z T-2 (s m-1)^n_sliding ~> Pa (s m-1)^n_sliding]
    allocate(CS%coef_prefactor(isd:ied,jsd:jed), source=0.0)
    allocate(CS%fB_elem(isd:ied,jsd:jed), source=0.0)
    allocate(CS%OD_av(isd:ied,jsd:jed), source=0.0)
    allocate(CS%ground_frac(isd:ied,jsd:jed), source=0.0)
    allocate(CS%taudx_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%taudy_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%sx_shelf(isd:ied,jsd:jed), source=0.0)
    allocate(CS%sy_shelf(isd:ied,jsd:jed), source=0.0)
    allocate(CS%bed_elev(isd:ied,jsd:jed), source=0.0)
    allocate(CS%bed_node(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%h_nodal(isd:ied,jsd:jed,1:2,1:2), source=0.0)
    allocate(CS%h_flot(isd:ied,jsd:jed,1:2,1:2), source=0.0)
    allocate(CS%Minv_xi(isd:ied,jsd:jed,1:2,1:2), source=0.0)
    allocate(CS%Minv_eta(isd:ied,jsd:jed,1:2,1:2), source=0.0)
    allocate(CS%cell_mean_w(isd:ied,jsd:jed,1:2,1:2), source=0.0)
    allocate(CS%h_source_rate(isd:ied,jsd:jed), source=0.0)
    allocate(CS%h_source_rate_last(isd:ied,jsd:jed), source=0.0)
    allocate(CS%phi_x_FV(IsdB:IedB,jsd:jed), source=1.0)
    allocate(CS%phi_y_FV(isd:ied,JsdB:JedB), source=1.0)
    allocate(CS%dg_art_visc_coef_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_art_visc_coef_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%dg_art_visc_nu_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_art_visc_nu_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%dg_art_visc_excess_frac_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_art_visc_excess_frac_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%dg_art_visc_allow_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_art_visc_allow_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%dg_art_visc_cell_scale(isd:ied,jsd:jed), source=1.0)
    allocate(CS%dg_slow_idle_face_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_slow_idle_face_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%mu_lim_xi(isd:ied,jsd:jed),    source=0.0)
    allocate(CS%mu_lim_eta(isd:ied,jsd:jed),   source=0.0)
    allocate(CS%mu_lim_cross(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_lim_phi_xi(isd:ied,jsd:jed),     source=1.0)
    allocate(CS%dg_lim_phi_eta(isd:ied,jsd:jed),    source=1.0)
    allocate(CS%dg_lim_phi_cross(isd:ied,jsd:jed),  source=1.0)
    allocate(CS%dg_lim_mass_drift(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_lim_phi(isd:ied,jsd:jed),        source=1.0)
    allocate(CS%dg_lim_pk_factor(isd:ied,jsd:jed),  source=0.0)
    allocate(CS%u_bdry_val(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%v_bdry_val(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%u_face_mask_bdry(IsdB:IedB,JsdB:JedB), source=-2.0)
    allocate(CS%v_face_mask_bdry(IsdB:iedB,JsdB:JedB), source=-2.0)
    allocate(CS%h_bdry_val(isd:ied,jsd:jed), source=0.0)

   ! additional restarts for ice shelf state
    call register_restart_field(CS%u_shelf, "u_shelf", .false., restart_CS, &
                                "ice sheet/shelf u-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%v_shelf, "v_shelf", .false., restart_CS, &
                                "ice sheet/shelf v-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%u_bdry_val, "u_bdry_val", .false., restart_CS, &
                                "ice sheet/shelf boundary u-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%v_bdry_val, "v_bdry_val", .false., restart_CS, &
                                "ice sheet/shelf boundary v-velocity", &
                                units="m s-1", conversion=US%L_T_to_m_s, hor_grid='Bu')
    call register_restart_field(CS%u_face_mask_bdry, "u_face_mask_bdry", .false., restart_CS, &
                                "ice sheet/shelf boundary u-mask", "nondim", hor_grid='Bu')
    call register_restart_field(CS%v_face_mask_bdry, "v_face_mask_bdry", .false., restart_CS, &
                                "ice sheet/shelf boundary v-mask", "nondim", hor_grid='Bu')

    call register_restart_field(CS%OD_av, "OD_av", .true., restart_CS, &
                                "Average open ocean depth in a cell", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%ground_frac, "ground_frac", .true., restart_CS, &
                                "fractional degree of grounding", "nondim")
    call register_restart_field(CS%C_basal_friction, "C_basal_friction", .true., restart_CS, &
                                "basal sliding coefficients", "Pa (s m-1)^n_sliding", conversion=US%RLZ_T2_to_Pa)
    call register_restart_field(CS%AGlen_visc, "AGlen_visc", .true., restart_CS, &
                                "ice-stiffness parameter", "Pa-3 s-1")
    call register_restart_field(CS%h_bdry_val, "h_bdry_val", .false., restart_CS, &
                                "ice thickness at the boundary", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%bed_elev, "bed elevation", .true., restart_CS, &
                                "bed elevation", "m", conversion=US%Z_to_m)
    ! Storage convention: h_shelf is the area-weighted cell mean Hbar derived from
    ! h_nodal; h_nodal(:,:,a,b) holds the 4 Q1 nodal corner values per cell.
    call register_restart_field(CS%h_nodal(:,:,1,1), "h_nodal_SW_DG", .true., restart_CS, &
                                "DG(1) nodal Q1 thickness, SW corner", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%h_nodal(:,:,2,1), "h_nodal_SE_DG", .true., restart_CS, &
                                "DG(1) nodal Q1 thickness, SE corner", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%h_nodal(:,:,1,2), "h_nodal_NW_DG", .true., restart_CS, &
                                "DG(1) nodal Q1 thickness, NW corner", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%h_nodal(:,:,2,2), "h_nodal_NE_DG", .true., restart_CS, &
                                "DG(1) nodal Q1 thickness, NE corner", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%first_dir_restart_IS, "first_direction_IS", .false., restart_CS, &
                                "Indicator of the first direction in split ice shelf calculations.", "nondim")
  endif

end subroutine register_ice_shelf_dyn_restarts

!> Initializes shelf model data, parameters and diagnostics
subroutine initialize_ice_shelf_dyn(param_file, Time, ISS, CS, G, US, diag, new_sim, Cp_ice, &
                                    Input_start_time, directory, solo_ice_sheet_in)
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(time_type),         intent(inout) :: Time !< The clock that that will indicate the model time
  type(ice_shelf_state),   intent(in)    :: ISS  !< A structure with elements that describe
                                                 !! the ice-shelf state
  type(ice_shelf_dyn_CS),  pointer       :: CS   !< A pointer to the ice shelf dynamics control structure
  type(ocean_grid_type),   intent(inout) :: G    !< The grid type describing the ice shelf grid.
  type(unit_scale_type),   intent(in)    :: US   !< A structure containing unit conversion factors
  type(diag_ctrl), target, intent(in)    :: diag !< A structure that is used to regulate the diagnostic output.
  logical,                 intent(in)    :: new_sim !< If true this is a new simulation, otherwise
                                                 !! has been started from a restart file.
  real,                    intent(in)    :: Cp_ice !< Heat capacity of ice [Q C-1 ~> J kg-1 degC-1]
  type(time_type),         intent(in)    :: Input_start_time !< The start time of the simulation.
  character(len=*),        intent(in)    :: directory  !< The directory where the ice sheet energy file goes.
  logical,       optional, intent(in)    :: solo_ice_sheet_in !< If present, this indicates whether
                                                 !! a solo ice-sheet driver.

  ! Local variables
  real    :: T_shelf_bdry ! A default ice shelf temperature to use for ice flowing
                          ! in through open boundaries [C ~> degC]
  !This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=200) :: IC_file,filename,inputdir
  character(len=40)  :: var_name
  character(len=40)  :: mdl = "MOM_ice_shelf_dyn"  ! This module's name.
  logical :: shelf_mass_is_dynamic, override_shelf_movement, active_shelf_dynamics
  logical :: enable_bugs  ! If true, the defaults for recently added bug-fix flags are set to
                          ! recreate the bugs, or if false bugs are only used if actively selected.
  logical :: debug
  integer :: i, j, isd, ied, jsd, jed, Isdq, Iedq, Jsdq, Jedq, iters
  logical :: valid_E, valid_W, valid_N, valid_S ! DG cold-start neighbour-mask checks
  real :: h_E, h_W, h_N, h_S ! Effective neighbour cell-mean thickness for DG cold-start [Z ~> m]
  logical :: slopes_from_file ! True if DG h_x, h_y were read from ICE_THICKNESS_FILE
  logical :: node_ic_used     ! True if nodal h was read and used to set h_shelf, h_x, h_y
  character(len=200) :: IS_energyfile  ! The name of the energy file.
  character(len=32) :: filename_appendix = '' ! FMS appendix to filename for ensemble runs
  character(len=16) :: inner_solver_str ! The type of inner solver to use for the SSA

  Isdq = G%isdB ; Iedq = G%iedB ; Jsdq = G%jsdB ; Jedq = G%jedB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) then
    call MOM_error(FATAL, "MOM_ice_shelf_dyn.F90, initialize_ice_shelf_dyn: "// &
                          "called with an associated control structure.")
    return
  endif
  if (CS%module_is_initialized) then
    call MOM_error(WARNING, "MOM_ice_shelf_dyn.F90, initialize_ice_shelf_dyn was "//&
             "called with a control structure that has already been initialized.")
  endif
  CS%module_is_initialized = .true.

  CS%diag => diag ! ; CS%Time => Time

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DEBUG", debug, default=.false.)
  call get_param(param_file, mdl, "DEBUG_IS", CS%debug, &
                 "If true, write verbose debugging messages for the ice shelf.", &
                 default=debug)
  call get_param(param_file, mdl, "DYNAMIC_SHELF_MASS", shelf_mass_is_dynamic, &
                 "If true, the ice sheet mass can evolve with time.", &
                 default=.false.)
  override_shelf_movement = .false. ; active_shelf_dynamics = .false.
  if (shelf_mass_is_dynamic) then
    call get_param(param_file, mdl, "OVERRIDE_SHELF_MOVEMENT", override_shelf_movement, &
                 "If true, user provided code specifies the ice-shelf "//&
                 "movement instead of the dynamic ice model.", default=.false., do_not_log=.true.)
    active_shelf_dynamics = .not.override_shelf_movement

    call get_param(param_file, mdl, "GROUNDING_LINE_INTERPOLATE", CS%GL_regularize, &
                 "If true, regularize the floatation condition at the "//&
                 "grounding line as in Goldberg Holland Schoof 2009.", default=.false.)
    call get_param(param_file, mdl, "GROUNDING_LINE_INTERP_SUBGRID_N", CS%n_sub_regularize, &
                 "The number of sub-partitions of each cell over which to "//&
                 "integrate for the interpolated grounding line. Each cell "//&
                 "is divided into NxN equally-sized rectangles, over which the "//&
                 "basal contribution is integrated by iterative quadrature.", &
                 default=0)
    call get_param(param_file, mdl, "GROUNDING_LINE_COUPLE", CS%GL_couple, &
                 "If true, let the floatation condition be determined by "//&
                 "ocean column thickness. This means that update_OD_ffrac "//&
                 "will be called.  GL_REGULARIZE and GL_COUPLE are exclusive.", &
                 default=.false., do_not_log=CS%GL_regularize)
    if (CS%GL_regularize) CS%GL_couple = .false.
    if (present(solo_ice_sheet_in)) then
      if (solo_ice_sheet_in) CS%GL_couple = .false.
    endif
    if (CS%GL_regularize .and. (CS%n_sub_regularize == 0)) call MOM_error (FATAL, &
      "GROUNDING_LINE_INTERP_SUBGRID_N must be a positive integer if GL regularization is used")
    call get_param(param_file, mdl, "ICE_SHELF_CFL_FACTOR", CS%CFL_factor, &
                 "A factor used to limit timestep as CFL_FACTOR * min (\Delta x / u). "//&
                 "This is only used with an ice-only model.", units="nondim", default=0.25)
  endif
  call get_param(param_file, mdl, "RHO_0", CS%density_ocean_avg, &
                 "avg ocean density used in floatation cond", &
                 units="kg m-3", default=1035., scale=US%kg_m3_to_R)
  if (active_shelf_dynamics) then
    call get_param(param_file, mdl, "ICE_VELOCITY_TIMESTEP", CS%velocity_update_time_step, &
                 "seconds between ice velocity calcs", units="s", scale=US%s_to_T, &
                 fail_if_missing=.true.)
    call get_param(param_file, mdl, "G_EARTH", CS%g_Earth, &
                 "The gravitational acceleration of the Earth.", &
                 units="m s-2", default=9.80, scale=US%m_s_to_L_T**2*US%Z_to_m)

    call get_param(param_file, mdl, "MIN_H_SHELF", CS%min_h_shelf, &
                 "min. ice thickness used during ice dynamics", &
                  units="m", default=0.,scale=US%m_to_Z)
    call get_param(param_file, mdl, "MIN_BASAL_TRACTION", CS%min_basal_traction, &
                 "min. allowed basal traction. Input is in [Pa m-1 yr], but is converted when read in to [Pa m-1 s]", &
                 units="Pa m-1 yr", default=0., scale=365.0*86400.0*US%Pa_to_RLZ_T2*US%L_T_to_m_s)
    call get_param(param_file, mdl, "MAX_SURFACE_SLOPE", CS%max_surface_slope, &
                 "max. allowed ice-sheet surface slope. To ignore, set to zero.", &
                 units="none", default=0., scale=US%m_to_Z/US%m_to_L)
    call get_param(param_file, mdl, "MIN_ICE_VISC", CS%min_ice_visc, &
                 "min. allowed Glen's law ice viscosity", &
                 units="Pa s", default=0., scale=US%Pa_to_RL2_T2*US%s_to_T)

    call get_param(param_file, mdl, "GLEN_EXPONENT", CS%n_glen, &
                 "nonlinearity exponent in Glen's Law", &
                  units="none", default=3.)
    call get_param(param_file, mdl, "MIN_STRAIN_RATE_GLEN", CS%eps_glen_min, &
                 "min. strain rate to avoid infinite Glen's law viscosity", &
                 units="s-1", default=1.e-19, scale=US%T_to_s)
    call get_param(param_file, mdl, "BASAL_FRICTION_EXP", CS%n_basal_fric, &
                 "Exponent in sliding law \tau_b = C u^(n_basal_fric)", &
                 units="none", fail_if_missing=.true.)
    call get_param(param_file, mdl, "USE_COULOMB_FRICTION", CS%CoulombFriction, &
                 "Use Coulomb Friction Law", &
                 units="none", default=.false., fail_if_missing=.false.)
    call get_param(param_file, mdl, "CF_MinN", CS%CF_MinN, &
                 "Minimum Coulomb friction effective pressure", &
                 units="Pa", default=1.0, scale=US%Pa_to_RLZ_T2, fail_if_missing=.false.)
    call get_param(param_file, mdl, "CF_PostPeak", CS%CF_PostPeak, &
                 "Coulomb friction post peak exponent", &
                 units="none", default=1.0, fail_if_missing=.false.)
    call get_param(param_file, mdl, "CF_Max", CS%CF_Max, &
                 "Coulomb friction maximum coefficient", &
                 units="none", default=0.5, fail_if_missing=.false.)
    ! Pre-compute Coulomb prefactor alpha = (q-1)^(q-1)/q^q for q=CF_PostPeak [nondim].
    ! Default is 1.0; only update when Coulomb is active and q /= 1.
    if (CS%CoulombFriction .and. CS%CF_PostPeak /= 1.0) &
      CS%alpha_coulomb = (CS%CF_PostPeak-1.0)**(CS%CF_PostPeak-1.0) / CS%CF_PostPeak**CS%CF_PostPeak

    call get_param(param_file, mdl, "DENSITY_ICE", CS%density_ice, &
                 "A typical density of ice.", units="kg m-3", default=917.0, scale=US%kg_m3_to_R)
    call get_param(param_file, mdl, "CONJUGATE_GRADIENT_TOLERANCE", CS%cg_tolerance, &
                 "For Picard iterations, the tolerance in CG solver, relative to initial residual", &
                 units="nondim", default=1.e-6)
    call get_param(param_file, mdl, "NEWTON_CONJUGATE_GRADIENT_TOLERANCE", CS%cg_newton_tolerance, &
                 "For inexact Newton iterations, the initial tolerance in CG solver, relative to initial residual", &
                 units="nondim", default=CS%cg_tolerance)
    CS%cg_tol_current = CS%cg_tolerance  ! Can be tightened adaptively during inexact Newton iterations
    call get_param(param_file, mdl, "ICE_NONLINEAR_TOLERANCE", CS%nonlinear_tolerance, &
                "nonlin tolerance in iterative velocity solve", units="nondim", default=1.e-6)
    call get_param(param_file, mdl, "NEWTON_AFTER_TOLERANCE", CS%newton_after_tolerance, &
                "Switch from Picard to Newton iterations in the nonlinear ice velocity solve when "//&
                "the fractional nonlinear residual falls below this tolerance. If <=0, no Picard.",&
                units="none", default=CS%nonlinear_tolerance)
    call get_param(param_file, mdl, "NEWTON_ADAPT_CG_TOL", CS%newton_adapt_cg_tol, &
                "Use an adaptive CG tolerance during Newton iterations.", default=.true.)
    call get_param(param_file, mdl, "NEWTON_EW_GAMMA", CS%ew_gamma, &
                "Gamma in Eisenstat-Walker adaptive Newton tolerance", units="nondim", default=0.9, &
                do_not_log=(.not. CS%newton_adapt_cg_tol))
    call get_param(param_file, mdl, "NEWTON_EW_ALPHA", CS%ew_alpha, &
                "Alpha in Eisenstat-Walker adaptive Newton tolerance", units="nondim", default=2.0,  &
                do_not_log=(.not. CS%newton_adapt_cg_tol))
    call get_param(param_file, mdl, "NEWTON_EW_SAFETY", CS%ew_safety, &
                "Safeguard Eisenstat-Walker using (0) no safeguard, (1) EW choice 2 threshold "//&
                "or (2) PETSc option 3 (Chacon 2008)", default=2, do_not_log=(.not. CS%newton_adapt_cg_tol))
    call get_param(param_file, mdl, "NEWTON_EW_1_THRESHOLD", CS%ew_1_thres, &
                "Eisenstat-Walker version 1 threshold", &
                units="nondim", default=0.1, do_not_log=(.not. CS%newton_adapt_cg_tol))
    call get_param(param_file, mdl, "NEWTON_EW_ETA_MAX", CS%ew_eta_max, &
                "Maximum allowed Eisenstat-Walker eta (between 0 and 1)", &
                units="nondim", default=0.9, do_not_log=(.not. CS%newton_adapt_cg_tol))
    if (CS%ew_eta_max<=0 .or. CS%ew_eta_max>= 1) &
      call MOM_error(FATAL, "NEWTON_EW_ETA_MAX must be between 0 and 1.")
    call get_param(param_file, mdl, "ICE_SHELF_INNER_SOLVER", inner_solver_str, &
                "Choice of inner linear solver for the ice-shelf SSA velocity system. "//&
                "Valid choices are CG (default), CR, and MINRES.", &
                default="CG")
    select case (trim(inner_solver_str))
      case ("CG")
        CS%inner_solver = INNER_CG
      case ("MINRES")
        CS%inner_solver = INNER_MINRES
      case ("CR")
        CS%inner_solver = INNER_CR
    end select
    call get_param(param_file, mdl, "CG_HALO_SHRINK", CS%cg_halo_shrink, &
                "If true, CG uses halo-shrinking to defer pass_vector calls. "//&
                "If false, uses a fixed CG_action range with one pass_vector(D) per iteration, "//&
                "which may reduce total communication for typical halo widths.", &
                default=.true.)
    call get_param(param_file, mdl, "CONJUGATE_GRADIENT_MAXIT", CS%cg_max_iterations, &
                "max iteratiions in CG solver", default=2000)
    call get_param(param_file, mdl, "THRESH_FLOAT_COL_DEPTH", CS%thresh_float_col_depth, &
                "min ocean thickness to consider ice *floating*; "//&
                "will only be important with use of tides", &
                units="m", default=1.e-3, scale=US%m_to_Z)
    call get_param(param_file, mdl, "NONLIN_SOLVE_ERR_MODE", CS%nonlin_solve_err_mode, &
                "Choose whether nonlin error in vel solve is based on nonlinear "//&
                "Linf norm residual (1), Linf norm relative change since last iteration (2), "//&
                "change in solution L2 norm (3), L2 norm residual (4), L2 backward norm (5)", default=3)
    if (CS%nonlin_solve_err_mode /= 5) then
      call get_param(param_file, mdl, "SSA_ADD_REL_RESID", CS%ssa_add_rel_resid, &
                  "Nonlinear error in vel solve will also depend on "// &
                  "L2 residual norm relative to RHS norm.", default=.false.)
    else
      CS%ssa_add_rel_resid = .false. !Avoids redundantly calculating err_mode 5 twice
    endif
    call get_param(param_file, mdl, "ICE_RR_NONLINEAR_TOLERANCE", CS%rr_nonlinear_tolerance, &
              "if ssa_add_rel_resid, the additional nonlin tolerance "//&
              "in the iterative velocity solve for the residual norm relative to RHS norm", &
              units="nondim", default=1.e-4)
    call get_param(param_file, mdl, "SHELF_MOVING_FRONT", CS%moving_shelf_front, &
                 "Specify whether to advance shelf front (and calve).", &
                 default=.false.)
    call get_param(param_file, mdl, "CALVE_TO_MASK", CS%calve_to_mask, &
                 "If true, do not allow an ice shelf where prohibited by a mask.", &
                 default=.false.)
    call get_param(param_file, mdl, "ADVECT_SHELF", CS%advect_shelf, &
                 "If true, advect ice shelf and evolve thickness", &
                 default=.true.)
    call get_param(param_file, mdl, "USE_DG_THICKNESS", CS%use_DG_thickness, &
                 "If true, use a DG(1) polynomial representation for ice thickness "//&
                 "with unsplit RK2 advection and sub-element Gauss quadrature for "//&
                 "driving stress. Requires h_x and h_y slope moments.", &
                 default=.false.)
    call get_param(param_file, mdl, "USE_NODAL_BED_FILE", CS%use_nodal_bed_file, &
                 "If true, read bed elevation directly at B-grid nodes from "//&
                 "NODAL_BED_FILE into CS%bed_node and derive the cell-centered "//&
                 "CS%bed_elev by bilinear averaging of the four surrounding "//&
                 "nodes. Skips reconstruct_bed_to_nodes and skips the "//&
                 "BED_TOPO_FILE read in initialize_ice_flow_from_file. "//&
                 "Requires USE_DG_THICKNESS=True.", &
                 default=.false., do_not_log=.not.CS%use_DG_thickness)
    if (CS%use_nodal_bed_file .and. .not. CS%use_DG_thickness) &
      call MOM_error(FATAL, "MOM_ice_shelf_dynamics: USE_NODAL_BED_FILE=True requires USE_DG_THICKNESS=True")
    call read_nodal_limiter_params(param_file, mdl, CS, US)
    call get_param(param_file, mdl, "REENTRANT_X", CS%reentrant_x, &
                 " If true, the domain is zonally reentrant.", &
                 default=.false.)
    call get_param(param_file, mdl, "REENTRANT_Y", CS%reentrant_y, &
                 " If true, the domain is meridionally reentrant.", &
                 default=.false.)
    call get_param(param_file, mdl, "ICE_VISCOSITY_COMPUTE", CS%ice_viscosity_compute, &
                 "If MODEL, compute ice viscosity internally using 1 or 4 quadrature points, "//&
                 "if OBS read from a file, "//&
                 "if CONSTANT a constant value (for debugging).", &
                 default="MODEL")

    call get_param(param_file, mdl, "ENABLE_BUGS_BY_DEFAULT", enable_bugs, &
                 default=.true., do_not_log=.true.)  ! This is logged from MOM.F90.
    call get_param(param_file, mdl, "ICE_SHELF_TOP_SLOPE_BUG", CS%shelf_top_slope_bugs, &
                 "If true, use directionally inconsistent estimates of the grid spacing when "//&
                 "calculating the ice shelf surface slope, and underestimate slopes near the "//&
                 "edge of the ice shelf by a factor of 2.", default=enable_bugs)

    if ((CS%visc_qps/=1) .and. (trim(CS%ice_viscosity_compute) /= "MODEL")) then
      call MOM_error(FATAL, "NUMBER_OF_ICE_VISCOSITY_QUADRATURE_POINTS must be 1 unless ICE_VISCOSITY_COMPUTE==MODEL.")
    endif
    call get_param(param_file, mdl, "INFLOW_SHELF_TEMPERATURE", T_shelf_bdry, &
                 "A default ice shelf temperature to use for ice flowing in through "//&
                 "open boundaries.", units="degC", default=-15.0, scale=US%degC_to_C)
  endif
  call get_param(param_file, mdl, "MISSING_SHELF_TEMPERATURE", CS%T_shelf_missing, &
                 "An ice shelf temperature to use where there is no ice shelf.",&
                 units="degC", default=-10.0, scale=US%degC_to_C)
  call get_param(param_file, mdl, "MIN_THICKNESS_SIMPLE_CALVE", CS%min_thickness_simple_calve, &
                 "Min thickness rule for the VERY simple calving law",&
                 units="m", default=0.0, scale=US%m_to_Z)
  CS%Cp_ice = Cp_ice !Heat capacity of ice (J kg-1 K-1), needed for heat flux of any bergs calved from
                     !the ice shelf and for ice sheet temperature solver
  !for write_ice_shelf_energy
      ! Note that the units of CS%Timeunit are the MKS units of [s].
  call get_param(param_file, mdl, "TIMEUNIT", CS%Timeunit, &
    "The time unit in seconds a number of input fields", &
    units="s", default=86400.0)
  if (CS%Timeunit < 0.0) CS%Timeunit = 86400.0
  call get_param(param_file, mdl, "ENERGYSAVEDAYS",CS%energysavedays, &
    "The interval in units of TIMEUNIT between saves of the "//&
    "energies of the run and other globally summed diagnostics.",&
    default=set_time(0,days=1), timeunit=CS%Timeunit)
  call get_param(param_file, mdl, "ENERGYSAVEDAYS_GEOMETRIC",CS%energysavedays_geometric, &
    "The starting interval in units of TIMEUNIT for the first call "//&
    "to save the energies of the run and other globally summed diagnostics. "//&
    "The interval increases by a factor of 2. after each call to write_ice_shelf_energy.",&
    default=set_time(seconds=0), timeunit=CS%Timeunit)
  if ((time_type_to_real(CS%energysavedays_geometric) > 0.) .and. &
    (CS%energysavedays_geometric < CS%energysavedays)) then
    CS%energysave_geometric = .true.
  else
    CS%energysave_geometric = .false.
  endif
  CS%Start_time = Input_start_time
  call get_param(param_file, mdl, "ICE_SHELF_ENERGYFILE", IS_energyfile, &
                 "The file to use to write the energies and globally "//&
                 "summed diagnostics.", default="ice_shelf.stats")
  !query fms_io if there is a filename_appendix (for ensemble runs)
  call get_filename_appendix(filename_appendix)
  if (len_trim(filename_appendix) > 0) then
    IS_energyfile = trim(IS_energyfile) //'.'//trim(filename_appendix)
  endif

  CS%IS_energyfile = trim(slasher(directory))//trim(IS_energyfile)
  call log_param(param_file, mdl, "output_path/ENERGYFILE", CS%IS_energyfile)
#ifdef STATSLABEL
  CS%IS_energyfile = trim(CS%IS_energyfile)//"."//trim(adjustl(STATSLABEL))
#endif

  ! Allocate memory in the ice shelf dynamics control structure that was not
  ! previously allocated for registration for restarts.

  if (active_shelf_dynamics) then
    allocate( CS%t_bdry_val(isd:ied,jsd:jed), source=T_shelf_bdry) ! [C ~> degC]
    allocate( CS%u_face_mask(Isdq:Iedq,Jsdq:Jedq), source=0.0)
    allocate( CS%v_face_mask(Isdq:Iedq,Jsdq:Jedq), source=0.0)
    allocate( CS%u_flux_bdry_val(Isdq:Iedq,jsd:jed), source=0.0)
    allocate( CS%v_flux_bdry_val(isd:ied,Jsdq:Jedq), source=0.0)
    allocate( CS%umask(Isdq:Iedq,Jsdq:Jedq), source=-1.0)
    allocate( CS%vmask(Isdq:Iedq,Jsdq:Jedq), source=-1.0)
    allocate( CS%tmask(Isdq:Iedq,Jsdq:Jedq), source=-1.0)

    CS%OD_rt_counter = 0
    allocate( CS%OD_rt(isd:ied,jsd:jed), source=0.0)
    allocate( CS%ground_frac_rt(isd:ied,jsd:jed), source=0.0)

    if (CS%calve_to_mask) then
      allocate( CS%calve_mask(isd:ied,jsd:jed), source=0.0)
    endif

    allocate(CS%Phi(1:8,1:4,isd:ied,jsd:jed), source=0.0)
    allocate(CS%Jac(1:4,isd:ied,jsd:jed), source=0.0)
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      call bilinear_shape_fn_grid(G, i, j, CS%Phi(:,:,i,j), CS%Jac(:,i,j))
    enddo ; enddo

    ! Per-cell nodal DG(1) metric tables (Minv_xi, Minv_eta, cell_mean_w).
    call init_nodal_DG_metric(CS, G)

    if (CS%GL_regularize) then
      allocate(CS%Phisub(2,2,CS%n_sub_regularize,CS%n_sub_regularize,2,2), source=0.0)
      call bilinear_shape_functions_subgrid(CS%Phisub, CS%n_sub_regularize)
    endif

    if ((trim(CS%ice_viscosity_compute) == "MODEL") .and. CS%visc_qps==1) then
      !for calculating viscosity and 1 cell-centered quadrature point per cell
      allocate(CS%PhiC(1:8,G%isc:G%iec,G%jsc:G%jec), source=0.0)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        call bilinear_shape_fn_grid_1qp(G, i, j, CS%PhiC(:,i,j))
      enddo ; enddo
    endif

    CS%elapsed_velocity_time = 0.0

    call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)
  endif

  ! Take additional initialization steps, for example of dependent variables.
  if (active_shelf_dynamics .and. .not.new_sim) then

    call pass_var(CS%OD_av,G%domain, complete=.false.)
    call pass_var(CS%ground_frac, G%domain, complete=.false.)
    call pass_var(CS%AGlen_visc, G%domain, complete=.false.)
    call pass_var(CS%bed_elev, G%domain, complete=.false.)
    call pass_var(CS%C_basal_friction, G%domain, complete=.false.)
    call pass_var(CS%h_bdry_val, G%domain, complete=.true.)
    call pass_var(CS%ice_visc, G%domain)
    if (CS%use_DG_thickness) then
      if (CS%use_nodal_bed_file) then
        call initialize_bed_node_from_file(CS%bed_node, CS%bed_elev, G, US, param_file)
      else
        call reconstruct_bed_to_nodes(CS, G, ISS%hmask)
      endif
    endif

    call pass_vector(CS%u_bdry_val, CS%v_bdry_val, G%domain, TO_ALL, BGRID_NE, complete=.false.)
    call pass_vector(CS%u_face_mask_bdry, CS%v_face_mask_bdry, G%domain, TO_ALL, BGRID_NE, complete=.true.)
    call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)

    ! This is unfortunately necessary (?); if grid is not symmetric the boundary values
    ! of u and v are otherwise not set till the end of the first linear solve, and so
    ! viscosity is not calculated correctly.
    ! This has to occur after init_boundary_values or some of the arrays on the
    ! right hand side have not been set up yet.
    if (.not. G%symmetric) then
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        if ((i+G%idg_offset) == (G%domain%nihalo+1)) then
          if (CS%u_face_mask(I-1,j) == 3) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I-1,J) = CS%u_bdry_val(I-1,J)
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I-1,J) = CS%v_bdry_val(I-1,J)
          elseif (CS%u_face_mask(I-1,j) == 5) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I-1,J) = CS%u_bdry_val(I-1,J)
          elseif (CS%u_face_mask(I-1,j) == 6) then
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I-1,J) = CS%v_bdry_val(I-1,J)
          endif
        endif
        if ((j+G%jdg_offset) == (G%domain%njhalo+1)) then
          if (CS%v_face_mask(i,J-1) == 3) then
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I,J-1) = CS%v_bdry_val(I,J-1)
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I,J-1) = CS%u_bdry_val(I,J-1)
          elseif (CS%v_face_mask(i,J-1) == 5) then
            CS%v_shelf(I-1,J-1) = CS%v_bdry_val(I-1,J-1)
            CS%v_shelf(I,J-1) = CS%v_bdry_val(I,J-1)
          elseif (CS%v_face_mask(i,J-1) == 6) then
            CS%u_shelf(I-1,J-1) = CS%u_bdry_val(I-1,J-1)
            CS%u_shelf(I,J-1) = CS%u_bdry_val(I,J-1)
          endif
        endif
      enddo ; enddo
    endif
    call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE)
  endif

  if (active_shelf_dynamics) then
    if (CS%first_dir_restart_IS > -1.0) then
      CS%first_direction_IS = modulo(NINT(CS%first_dir_restart_IS), 2)
    else
      CS%first_dir_restart_IS = real(modulo(CS%first_direction_IS, 2))
    endif

    ! If we are calving to a mask, i.e. if a mask exists where a shelf cannot, read the mask from a file.
    if (CS%calve_to_mask) then
      call MOM_mesg("  MOM_ice_shelf.F90, initialize_ice_shelf: reading calving_mask")

      call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
      inputdir = slasher(inputdir)
      call get_param(param_file, mdl, "CALVING_MASK_FILE", IC_file, &
                   "The file with a mask for where calving might occur.", &
                   default="ice_shelf_h.nc")
      call get_param(param_file, mdl, "CALVING_MASK_VARNAME", var_name, &
                   "The variable to use in masking calving.", &
                   default="area_shelf_h")

      filename = trim(inputdir)//trim(IC_file)
      call log_param(param_file, mdl, "INPUTDIR/CALVING_MASK_FILE", filename)
      if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
         " calving mask file: Unable to open "//trim(filename))

      call MOM_read_data(filename,trim(var_name),CS%calve_mask,G%Domain)
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        if (CS%calve_mask(i,j) > 0.0) CS%calve_mask(i,j) = 1.0
      enddo ; enddo
      call pass_var(CS%calve_mask,G%domain)
    endif

    ! initialize basal friction coefficients
    if (new_sim) then
      call initialize_ice_C_basal_friction(CS%C_basal_friction, G, US, param_file)
      call pass_var(CS%C_basal_friction, G%domain, complete=.false.)

      ! initialize ice-stiffness AGlen
      call initialize_ice_AGlen(CS%AGlen_visc, CS%ice_viscosity_compute, G, US, param_file)
      call pass_var(CS%AGlen_visc, G%domain, complete=.false.)

      !initialize boundary conditions
      call initialize_ice_shelf_boundary_from_file(CS%u_face_mask_bdry, CS%v_face_mask_bdry, &
                  CS%u_bdry_val, CS%v_bdry_val, CS%umask, CS%vmask, CS%h_bdry_val, &
                  ISS%hmask,  ISS%h_shelf, G, US, param_file )
      call pass_var(ISS%hmask, G%domain, complete=.false.)
      call pass_var(CS%h_bdry_val, G%domain, complete=.true.)
      call pass_vector(CS%u_bdry_val, CS%v_bdry_val, G%domain, TO_ALL, BGRID_NE, complete=.false.)
      call pass_vector(CS%u_face_mask_bdry, CS%v_face_mask_bdry, G%domain, TO_ALL, BGRID_NE, complete=.false.)

      !initialize ice flow characteristic (velocities, bed elevation under the grounded part, etc) from file
      call initialize_ice_flow_from_file(CS%bed_elev,CS%u_shelf, CS%v_shelf, CS%ground_frac, &
                  G, US, param_file, skip_bed=CS%use_nodal_bed_file)
      call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE, complete=.true.)
      call pass_var(CS%ground_frac, G%domain, complete=.true.)
      if (CS%use_nodal_bed_file) then
        ! Reads bed_node from file and derives bed_elev (both halo-updated inside).
        call initialize_bed_node_from_file(CS%bed_node, CS%bed_elev, G, US, param_file)
      else
        call pass_var(CS%bed_elev, G%domain, complete=.true.)
        if (CS%use_DG_thickness) call reconstruct_bed_to_nodes(CS, G, ISS%hmask)
      endif
      if (CS%use_DG_thickness) then
        ! Nodal DG(1) cold-start: try node-file IC, otherwise project cell means.
        call pass_var(ISS%h_shelf, G%domain)
        CS%h_nodal(:,:,:,:) = 0.0
        node_ic_used = .false.
        call initialize_DG_thickness_from_node_file(ISS%h_shelf, CS%h_nodal, ISS%hmask, &
                                                    node_ic_used, G, US, param_file)
        if (.not. node_ic_used) then
          call initialize_h_nodal_from_cellmean(ISS%h_shelf, CS%h_nodal, ISS%hmask, G)
        endif
        if (CS%nodal_positivity) call nodal_positivity_limit(CS, G, ISS)
        call pass_corner_field(CS%h_nodal, G)
        call enforce_wrap_corner_consistency(CS, ISS, G)
        call recompute_h_shelf_from_nodal(CS, ISS, G)
      endif
      call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)

      do J=Jsdq,Jedq ; do I=Isdq,Iedq
        if (CS%umask(I,J) == 3) then
          CS%u_shelf(I,J) = CS%u_bdry_val(I,J)
        elseif (CS%umask(I,J) == 0) then
          CS%u_shelf(I,J) = 0
        endif
        if (CS%vmask(I,J) == 3) then
          CS%v_shelf(I,J) = CS%v_bdry_val(I,J)
        elseif (CS%vmask(I,J) == 0) then
          CS%v_shelf(I,J) = 0
        endif
      enddo ; enddo
    endif

  ! Register diagnostics.
    CS%id_u_shelf = register_diag_field('ice_shelf_model','u_shelf',CS%diag%axesB1, Time, &
       'x-velocity of ice', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_v_shelf = register_diag_field('ice_shelf_model','v_shelf',CS%diag%axesB1, Time, &
       'y-velocity of ice', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_shelf_speed = register_diag_field('ice_shelf_model','shelf_speed',CS%diag%axesB1, Time, &
       'speed of of ice shelf', 'm yr-1', conversion=365.0*86400.0*US%L_T_to_m_s)
    CS%id_taudx_shelf = register_diag_field('ice_shelf_model','taudx_shelf',CS%diag%axesB1, Time, &
       'x-driving stress of ice', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_taudy_shelf = register_diag_field('ice_shelf_model','taudy_shelf',CS%diag%axesB1, Time, &
       'y-driving stress of ice', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_taud_shelf = register_diag_field('ice_shelf_model','taud_shelf',CS%diag%axesB1, Time, &
       'magnitude of driving stress of ice', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_sx_shelf = register_diag_field('ice_shelf_model', 'sx_shelf', CS%diag%axesT1, Time, &
       'x-surface slope of ice', 'none')
    CS%id_sy_shelf = register_diag_field('ice_shelf_model', 'sy_shelf', CS%diag%axesT1, Time, &
       'y-surface slope of ice', 'none')
    CS%id_surf_slope_mag_shelf = register_diag_field('ice_shelf_model', 'surf_slope_mag_shelf', CS%diag%axesT1, Time, &
       'magnitude of surface slope of ice', 'none')
    CS%id_u_mask = register_diag_field('ice_shelf_model','u_mask',CS%diag%axesB1, Time, &
       'mask for u-nodes', 'none')
    CS%id_v_mask = register_diag_field('ice_shelf_model','v_mask',CS%diag%axesB1, Time, &
       'mask for v-nodes', 'none')
    CS%id_ground_frac = register_diag_field('ice_shelf_model','ice_ground_frac',CS%diag%axesT1, Time, &
       'fraction of cell that is grounded; under GL_regularize this is the fraction of '//&
       'sub-cell quadrature points whose draft sits below the bed', 'none')
    CS%id_col_thick = register_diag_field('ice_shelf_model','col_thick',CS%diag%axesT1, Time, &
       'ocean column thickness passed to ice model', 'm', conversion=US%Z_to_m)
    CS%id_visc_shelf = register_diag_field('ice_shelf_model','ice_visc',CS%diag%axesT1, Time, &
       'vi-viscosity', 'Pa m s', conversion=US%RL2_T2_to_Pa*US%Z_to_m*US%T_to_s) !vertically integrated viscosity
    CS%id_taub = register_diag_field('ice_shelf_model','taub_beta',CS%diag%axesT1, Time, &
       'taub', units='MPa yr m-1', conversion=1e-6*US%RLZ_T2_to_Pa/(365.0*86400.0*US%L_T_to_m_s))
    CS%id_OD_av = register_diag_field('ice_shelf_model','OD_av',CS%diag%axesT1, Time, &
       'intermediate ocean column thickness passed to ice model', 'm', conversion=US%Z_to_m)

    if (CS%use_DG_thickness) then
      CS%id_bed_node = register_diag_field('ice_shelf_model','bed_node',CS%diag%axesB1, Time, &
         'Bed elevation at B-grid nodes (DG bilinear reconstruction)', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_SW = register_diag_field('ice_shelf_model','h_nodal_SW',CS%diag%axesT1, Time, &
         'DG(1) nodal Q1 thickness at SW cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_SE = register_diag_field('ice_shelf_model','h_nodal_SE',CS%diag%axesT1, Time, &
         'DG(1) nodal Q1 thickness at SE cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_NW = register_diag_field('ice_shelf_model','h_nodal_NW',CS%diag%axesT1, Time, &
         'DG(1) nodal Q1 thickness at NW cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_NE = register_diag_field('ice_shelf_model','h_nodal_NE',CS%diag%axesT1, Time, &
         'DG(1) nodal Q1 thickness at NE cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_jump_node = register_diag_field('ice_shelf_model','h_jump_node',CS%diag%axesB1, Time, &
         'DG(1) max-minus-min of co-located corner thickness across up to 4 touching cells '//&
         '(hmask=1 only) at each B-grid node', 'm', conversion=US%Z_to_m)
      CS%id_h_jump_node_rel = register_diag_field('ice_shelf_model','h_jump_node_rel',CS%diag%axesB1, Time, &
         'DG(1) B-node jump as a fraction of the mean cell-mean thickness over the touching cells '//&
         '(hmask=1 only)', 'nondim')
      CS%id_h_node_max = register_diag_field('ice_shelf_model','h_node_max',CS%diag%axesB1, Time, &
         'DG(1) maximum co-located corner thickness across up to 4 touching cells (hmask=1 only) '//&
         'at each B-grid node. Paired with h_node_min, separates a high-side spike (h_node_max '//&
         'large vs neighbor cell means) from a low-side pit (h_node_min small) at the same node.', &
         'm', conversion=US%Z_to_m)
      CS%id_h_node_min = register_diag_field('ice_shelf_model','h_node_min',CS%diag%axesB1, Time, &
         'DG(1) minimum co-located corner thickness across up to 4 touching cells (hmask=1 only) '//&
         'at each B-grid node. See h_node_max.', 'm', conversion=US%Z_to_m)
      CS%id_h_jump_envelope = register_diag_field('ice_shelf_model','h_jump_envelope', &
         CS%diag%axesB1, Time, &
         'Per-B-node cell-mean envelope width Hmax_B - Hmin_B, built over hmask=1 and hmask=3 '//&
         '(Dirichlet) cells touching the node. Reports the local roughness of the cell-mean '//&
         'thickness field around each B-node.', &
         'm', conversion=US%Z_to_m)
      CS%id_h_jump_envelope_rel = register_diag_field('ice_shelf_model','h_jump_envelope_rel', &
         CS%diag%axesB1, Time, &
         'h_jump_envelope normalised by the mean of the contributing cell means.', 'nondim')
      CS%id_h_overshoot_node = register_diag_field('ice_shelf_model','h_overshoot_node', &
         CS%diag%axesB1, Time, &
         'Per-B-node Barth-Jespersen overshoot: max over the touching DG(1) corner values of '//&
         'max(0, h_corner - Hmax_B, Hmin_B - h_corner), where [Hmin_B, Hmax_B] is the envelope '//&
         'of cell-mean thicknesses over the cells touching the node. Nonzero only where a DG '//&
         'corner has wandered outside the local neighbor-mean envelope (true sub-cell '//&
         'discontinuous mode, i.e. an unphysical overshoot rather than a faithful resolved '//&
         'sharp gradient).', &
         'm', conversion=US%Z_to_m)
      CS%id_h_overshoot_node_rel = register_diag_field('ice_shelf_model','h_overshoot_node_rel', &
         CS%diag%axesB1, Time, &
         'h_overshoot_node normalised by 0.5*(Hmax_B + Hmin_B).', 'nondim')
      CS%id_s_overshoot_node = register_diag_field('ice_shelf_model','s_overshoot_node', &
         CS%diag%axesB1, Time, &
         'Per-B-node Barth-Jespersen overshoot in the surface-elevation field: max over the '//&
         'touching DG(1) corner s values of max(0, s_corner - Smax_B, Smin_B - s_corner), with '//&
         'corner s computed by per-side flotation from h_corner and bed_node, and [Smin_B, '//&
         'Smax_B] the envelope of cell-mean surfaces (Hbar projected with cell-center bed via '//&
         'flotation) over the cells touching the node. The s-space analogue of h_overshoot_node: '//&
         'isolates the driving-stress-relevant spurious surface mode from the harmless '//&
         'thickness-only overshoot that disappears across the grounding line under flotation.', &
         'm', conversion=US%Z_to_m)
      CS%id_s_overshoot_node_rel = register_diag_field('ice_shelf_model','s_overshoot_node_rel', &
         CS%diag%axesB1, Time, &
         's_overshoot_node normalised by 0.5*(|Smax_B| + |Smin_B|).', 'nondim')
      CS%id_h_source_rate = register_diag_field('ice_shelf_model','h_source_rate',CS%diag%axesT1, Time, &
         'Cell-mean thickness source rate (basal melt + surface SMB) consumed by the last DG advect step', &
         'm s-1', conversion=US%Z_to_m*US%s_to_T)
      CS%id_dg_art_visc_coef_u = register_diag_field('ice_shelf_model','dg_art_visc_coef_u', &
         CS%diag%axesCu1, Time, &
         'Per-face DG(1) artificial-viscosity coefficient on u-faces (post per-cell CFL '//&
         'scaling) from the last spatial-operator call. Smooth faces report ~0 (smoothness '//&
         'gate off); shocky faces report up to DG1_ART_VISC_C_MAX; cells whose summed face '//&
         'rates would exceed the per-cell stability budget are scaled below c_max.', 'nondim')
      CS%id_dg_art_visc_coef_v = register_diag_field('ice_shelf_model','dg_art_visc_coef_v', &
         CS%diag%axesCv1, Time, &
         'Per-face DG(1) artificial-viscosity coefficient on v-faces. See dg_art_visc_coef_u.', &
         'nondim')
      CS%id_dg_art_visc_nu_u = register_diag_field('ice_shelf_model','dg_art_visc_nu_u', &
         CS%diag%axesCu1, Time, &
         'Per-face effective DG(1) artificial viscosity on u-faces (post per-cell CFL scaling), '//&
         'nu = c_face * u_eff_face_mean * dx_perp, with u_eff = |u_face| + '//&
         'DG1_ART_VISC_STRAIN_COEF * eps_e_face * dx_perp. Comparable to a physical '//&
         'diffusivity: face flux = nu * [h_eq] / dx_perp * dy_face. Distinguishes '//&
         'gate-active-but-quiet faces (low u_eff -> small nu despite c_face = c_max) from '//&
         'gate-active-and-damping faces (large u_eff -> large nu). Use with c_face to '//&
         'separate gate response from actual damping rate.', &
         'm2 s-1', conversion=US%L_T_to_m_s*US%L_to_m)
      CS%id_dg_art_visc_nu_v = register_diag_field('ice_shelf_model','dg_art_visc_nu_v', &
         CS%diag%axesCv1, Time, &
         'As dg_art_visc_nu_u but on v-faces.', &
         'm2 s-1', conversion=US%L_T_to_m_s*US%L_to_m)
      CS%id_dg_art_visc_excess_frac_u = register_diag_field('ice_shelf_model', &
         'dg_art_visc_excess_frac_u', CS%diag%axesCu1, Time, &
         'Fraction of the well-balanced surface jump treated as excess by the DG(1) '//&
         'artificial-viscosity EXCESS_JUMP allowance on u-faces (max over the 2 face QPs '//&
         'of |ds_use|/|ds_qp|). 1 with near-zero dg_art_visc_allow_u means the cell means '//&
         'lend the jump no support (slope-only structure, or a grounding-line-straddle '//&
         'face where the mean-branch allowance breaks down); 1 with a large exceeded '//&
         'allowance means genuine envelope escape (viscosity load-bearing); << 1 means '//&
         'only the edge of mean-supported structure is being shaved. Only meaningful '//&
         'with DG1_ART_VISC_EXCESS_JUMP=True (reports 0 otherwise).', 'nondim')
      CS%id_dg_art_visc_excess_frac_v = register_diag_field('ice_shelf_model', &
         'dg_art_visc_excess_frac_v', CS%diag%axesCv1, Time, &
         'As dg_art_visc_excess_frac_u but on v-faces.', 'nondim')
      CS%id_dg_art_visc_allow_u = register_diag_field('ice_shelf_model', &
         'dg_art_visc_allow_u', CS%diag%axesCu1, Time, &
         'Mean-supported surface-jump allowance |ds_bar| of the DG(1) artificial-'//&
         'viscosity EXCESS_JUMP band on u-faces (max over the 2 face QPs): the surface '//&
         'jump the two cell means imply at the face-QP bed (with '//&
         'DG1_ART_VISC_EXCESS_BRANCH_MAX, the max over admissible flotation-branch '//&
         'assignments of the means). Compare with s_jump_face_u '//&
         'and dg_art_visc_excess_frac_u to separate envelope damping from legitimate-'//&
         'structure shaving. Only meaningful with DG1_ART_VISC_EXCESS_JUMP=True '//&
         '(reports 0 otherwise).', 'm', conversion=US%Z_to_m)
      CS%id_dg_art_visc_allow_v = register_diag_field('ice_shelf_model', &
         'dg_art_visc_allow_v', CS%diag%axesCv1, Time, &
         'As dg_art_visc_allow_u but on v-faces.', 'm', conversion=US%Z_to_m)
      CS%id_dg_art_visc_cell_scale = register_diag_field('ice_shelf_model', &
         'dg_art_visc_cell_scale', CS%diag%axesT1, Time, &
         'Per-cell DG(1) artificial-viscosity cap throttle factor (1 = stability cap '//&
         'dormant; < 1 = cap engaged, all face coefficients of the cell rescaled by '//&
         'this value so the summed jump-mode decay rates stay within the SSP-RK2 '//&
         'budget DG1_ART_VISC_KCELL).', 'nondim')
      CS%id_dg_slow_idle_face_u = register_diag_field('ice_shelf_model','dg_slow_idle_face_u', &
         CS%diag%axesCu1, Time, &
         'DG(1) stagnant-jump indicator on u-faces (0 or 1). 1 where |u_face| and the SSA '//&
         'effective strain rate are both below internal tiny thresholds while the well-'//&
         'balanced equivalent jump exceeds a tolerance. Diagnostic for whether the '//&
         'velocity-magnitude + strain-rate damping channels both starve simultaneously; '//&
         'persistent non-zero values motivate adding a constant velocity floor.', 'nondim')
      CS%id_dg_slow_idle_face_v = register_diag_field('ice_shelf_model','dg_slow_idle_face_v', &
         CS%diag%axesCv1, Time, &
         'DG(1) stagnant-jump indicator on v-faces. See dg_slow_idle_face_u.', 'nondim')
      CS%id_h_jump_face_u = register_diag_field('ice_shelf_model','h_jump_face_u', &
         CS%diag%axesCu1, Time, &
         'DG(1) broken-Q1 thickness jump across u-faces (max |[h]| over the 2 face nodes). '//&
         'Localises the inter-element discontinuity to a single face, unlike the B-node '//&
         'max-min h_jump_node which conflates the up-to-4 faces meeting at a corner.', &
         'm', conversion=US%Z_to_m)
      CS%id_h_jump_face_v = register_diag_field('ice_shelf_model','h_jump_face_v', &
         CS%diag%axesCv1, Time, &
         'DG(1) broken-Q1 thickness jump across v-faces (max |[h]| over the 2 face nodes).', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_u = register_diag_field('ice_shelf_model','s_jump_face_u', &
         CS%diag%axesCu1, Time, &
         'DG(1) surface-elevation jump across u-faces (max |[s]| over the 2 face nodes, '//&
         'per-side flotation). This is the quantity the sub-grid driving stress consumes '//&
         '(jump_factor = rho*g*{h}*[s]); unlike [h] it accounts for the nonlinear '//&
         'grounded/floating thickness-to-surface map across the grounding line.', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_v = register_diag_field('ice_shelf_model','s_jump_face_v', &
         CS%diag%axesCv1, Time, &
         'DG(1) surface-elevation jump across v-faces (max |[s]| over the 2 face nodes, '//&
         'per-side flotation).', 'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_u_rel = register_diag_field('ice_shelf_model','s_jump_face_u_rel', &
         CS%diag%axesCu1, Time, &
         'Surface-elevation jump on u-faces relative to mean cell thickness, '//&
         '|[s]| / max(min_h_shelf, 0.5*(Hbar_A + Hbar_B)). Approximates the relative '//&
         'driving-stress contamination from spurious broken-Q1 jumps; healthy shelf '//&
         'regions should sit well below 0.05 (steep grounded slopes can tolerate more).', &
         'nondim')
      CS%id_s_jump_face_v_rel = register_diag_field('ice_shelf_model','s_jump_face_v_rel', &
         CS%diag%axesCv1, Time, &
         'As s_jump_face_u_rel but on v-faces.', 'nondim')
      CS%id_h_jump_face_u_signed = register_diag_field('ice_shelf_model','h_jump_face_u_signed', &
         CS%diag%axesCu1, Time, &
         'Signed DG(1) thickness jump on u-faces, mean over the 2 face nodes of '//&
         '(h_plus - h_minus) with plus = east cell, minus = west cell. Positive when the '//&
         'east-side cell is thicker. Pair adjacent signed jumps to detect 2dx oscillations '//&
         '(sign-flip pattern) vs. resolved gradients (consistent sign).', &
         'm', conversion=US%Z_to_m)
      CS%id_h_jump_face_v_signed = register_diag_field('ice_shelf_model','h_jump_face_v_signed', &
         CS%diag%axesCv1, Time, &
         'Signed DG(1) thickness jump on v-faces, mean over the 2 face nodes of '//&
         '(h_plus - h_minus) with plus = north cell, minus = south cell.', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_u_signed = register_diag_field('ice_shelf_model','s_jump_face_u_signed', &
         CS%diag%axesCu1, Time, &
         'Signed surface-elevation jump on u-faces, mean over the 2 face nodes of '//&
         '(s_plus - s_minus) with per-side flotation and plus = east cell. This is the '//&
         'quantity the DG(1) artificial viscosity actually responds to (via the '//&
         'well-balanced equivalent thickness jump); flips in sign cell-to-cell along a '//&
         'shear margin indicate under-damped oscillations.', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_v_signed = register_diag_field('ice_shelf_model','s_jump_face_v_signed', &
         CS%diag%axesCv1, Time, &
         'Signed surface-elevation jump on v-faces, mean over the 2 face nodes of '//&
         '(s_plus - s_minus) with per-side flotation and plus = north cell.', &
         'm', conversion=US%Z_to_m)
      CS%id_un_face_u = register_diag_field('ice_shelf_model','un_face_u', &
         CS%diag%axesCu1, Time, &
         'Face-normal ice speed |u.n| on u-faces (mean of the 2 endpoint B-node u_shelf '//&
         'values). Pair with h_jump_face_u / s_jump_face_u to test whether jumps accumulate '//&
         'at shear-margin faces where u.n ~ 0.', 'm s-1', conversion=US%L_T_to_m_s)
      CS%id_un_face_v = register_diag_field('ice_shelf_model','un_face_v', &
         CS%diag%axesCv1, Time, &
         'Face-normal ice speed |v.n| on v-faces (mean of the 2 endpoint B-node v_shelf '//&
         'values). Pair with h_jump_face_v / s_jump_face_v.', 'm s-1', conversion=US%L_T_to_m_s)
      CS%id_dg_eps_face_u = register_diag_field('ice_shelf_model','dg_eps_face_u', &
         CS%diag%axesCu1, Time, &
         'Effective SSA strain rate eps_e at u-face midpoints, computed by the same '//&
         'boundary-aware stencil used by the DG(1) artificial viscosity. Multiply by '//&
         'DG1_ART_VISC_STRAIN_COEF * dxCu to get the strain-driven velocity floor '//&
         '(u_floor = alpha*eps_e*dx_perp); compare with un_face_u to see where the floor '//&
         'beats advection. Posted unconditionally so it can be used to tune '//&
         'DG1_ART_VISC_STRAIN_COEF before turning the viscosity on.', &
         'yr-1', conversion=365.0*86400.0*US%s_to_T)
      CS%id_dg_eps_face_v = register_diag_field('ice_shelf_model','dg_eps_face_v', &
         CS%diag%axesCv1, Time, &
         'As dg_eps_face_u but on v-faces.', &
         'yr-1', conversion=365.0*86400.0*US%s_to_T)
      CS%id_dg_lim_phi_xi = register_diag_field('ice_shelf_model','dg_lim_phi_xi', &
         CS%diag%axesT1, Time, &
         'Per-cell xi-slope (east-west) mode scaling factor from the DG(1) hierarchical '//&
         'limiter [0,1]. 1 = no limiting, 0 = mode fully collapsed.', 'nondim')
      CS%id_dg_lim_phi_eta = register_diag_field('ice_shelf_model','dg_lim_phi_eta', &
         CS%diag%axesT1, Time, &
         'Per-cell eta-slope (north-south) mode scaling factor from the DG(1) hierarchical '//&
         'limiter [0,1].', 'nondim')
      CS%id_dg_lim_phi_cross = register_diag_field('ice_shelf_model','dg_lim_phi_cross', &
         CS%diag%axesT1, Time, &
         'Per-cell cross (saddle/twist) mode scaling factor from the DG(1) hierarchical '//&
         'limiter [0,1].', 'nondim')
      CS%id_dg_lim_mass_drift = register_diag_field('ice_shelf_model','dg_lim_mass_drift', &
         CS%diag%axesT1, Time, &
         'Per-cell change in cell-mean thickness produced by the DG(1) hierarchical limiter. '//&
         'Should be ~machine epsilon when the orthogonalised mode templates are correct.', &
         'm', conversion=US%Z_to_m)
      CS%id_dg_lim_phi = register_diag_field('ice_shelf_model','dg_lim_phi', &
         CS%diag%axesT1, Time, &
         'Per-cell DG(1) limiter strength [0,1]. 1 = no limiting at this cell, 0 = full '//&
         'collapse. For the isotropic single-phi variant this is the unique scaling factor; '//&
         'for the anisotropic per-mode variant this is the min over (phi_xi, phi_eta, phi_cross), '//&
         'i.e. the most-limiting direction.', 'nondim')
      CS%id_dg_lim_pk_factor = register_diag_field('ice_shelf_model','dg_lim_pk_factor', &
         CS%diag%axesT1, Time, &
         'Park-Kim MLP-u2 smooth-extrema indicator from the DG(1) hierarchical limiter. '//&
         '1 where the second-difference sign-consistency check across the 3-cell stencil in '//&
         'both x and y flags the cell as a smooth extremum (full slack granted); 0 where the '//&
         'check rejects and strict MLP-u2 vertex bounds apply.', 'nondim')
    else
      CS%id_phi_x_FV = register_diag_field('ice_shelf_model','phi_x_FV',CS%diag%axesCu1, Time, &
         'Van Leer slope-limiter factor at each u-face from ice_shelf_advect_thickness_x '//&
         '(1=no clip / smooth balanced, 0=full clip, 2=max compressive; faces where the limiter '//&
         'branch was not taken report 1.0)', 'nondim')
      CS%id_phi_y_FV = register_diag_field('ice_shelf_model','phi_y_FV',CS%diag%axesCv1, Time, &
         'Van Leer slope-limiter factor at each v-face from ice_shelf_advect_thickness_y '//&
         '(range [0,2])', 'nondim')
    endif

    CS%id_duHdx = register_diag_field('ice_shelf_model','duHdx',CS%diag%axesT1, Time, &
       'x-component of ice-sheet flux divergence', 'm yr-1', conversion=365.0*86400.0*US%Z_to_m*US%s_to_T)
    CS%id_dvHdy = register_diag_field('ice_shelf_model','dvHdy',CS%diag%axesT1, Time, &
       'y-component of ice-sheet flux divergence', 'm yr-1', conversion=365.0*86400.0*US%Z_to_m*US%s_to_T)
    CS%id_fluxdiv = register_diag_field('ice_shelf_model','fluxdiv',CS%diag%axesT1, Time, &
       'ice-sheet flux divergence', 'm yr-1', conversion=365.0*86400.0*US%Z_to_m*US%s_to_T)
    CS%id_strainrate_xx = register_diag_field('ice_shelf_model','strainrate_xx',CS%diag%axesT1, Time, &
       'x-component of ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_strainrate_yy = register_diag_field('ice_shelf_model','strainrate_yy',CS%diag%axesT1, Time, &
       'y-component of ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_strainrate_xy = register_diag_field('ice_shelf_model','strainrate_xy',CS%diag%axesT1, Time, &
       'xy-component of ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_pstrainrate_1 = register_diag_field('ice_shelf_model','pstrainrate_1',CS%diag%axesT1, Time, &
       'max principal horizontal ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_pstrainrate_2 = register_diag_field('ice_shelf_model','pstrainrate_2',CS%diag%axesT1, Time, &
       'min principal horizontal ice-shelf strain-rate', 'yr-1', conversion=365.0*86400.0*US%s_to_T)
    CS%id_devstress_xx = register_diag_field('ice_shelf_model','devstress_xx',CS%diag%axesT1, Time, &
       'x-component of ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_devstress_yy = register_diag_field('ice_shelf_model','devstress_yy',CS%diag%axesT1, Time, &
       'y-component of ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_devstress_xy = register_diag_field('ice_shelf_model','devstress_xy',CS%diag%axesT1, Time, &
       'xy-component of ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_pdevstress_1 = register_diag_field('ice_shelf_model','pdevstress_1',CS%diag%axesT1, Time, &
       'max principal horizontal ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)
    CS%id_pdevstress_2 = register_diag_field('ice_shelf_model','pdevstress_2',CS%diag%axesT1, Time, &
       'min principal ice-shelf deviatoric stress', 'kPa', conversion=1.e-3*US%RLZ_T2_to_Pa)

    !Update these variables so that they are nonzero in case
    !IS_dynamics_post_data is called before update_ice_shelf
    if (CS%id_taudx_shelf>0 .or. CS%id_taudy_shelf>0) then
      if (CS%use_DG_thickness) then
        if (CS%dg_gl_gate_continuous) call compute_h_flot(CS, ISS, G)
        if (CS%dg_driving_stress_IBP) then
          call calc_shelf_driving_stress_DG(CS, ISS, G, US, CS%taudx_shelf, CS%taudy_shelf, CS%OD_av)
        else
          call calc_shelf_driving_stress_DG_strong(CS, ISS, G, US, CS%taudx_shelf, CS%taudy_shelf, CS%OD_av)
        endif
      else
        call calc_shelf_driving_stress(CS, ISS, G, US, CS%taudx_shelf, CS%taudy_shelf, CS%OD_av)
      endif
    endif
    if (CS%id_visc_shelf>0) then
      call calc_shelf_visc(CS, ISS, G, US, CS%u_shelf, CS%v_shelf)
    endif
  endif

  if (new_sim) then
    call update_OD_ffrac_uncoupled(CS, G, ISS%h_shelf(:,:))
  endif

end subroutine initialize_ice_shelf_dyn


subroutine initialize_diagnostic_fields(CS, ISS, G, US, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G   !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US   !< A structure containing unit conversion factors
  type(time_type),        intent(in)    :: Time !< The current model time

  integer         :: i, j, iters, isd, ied, jsd, jed
  real            :: rhoi_rhow
  real            :: OD  ! Depth of open water below the ice shelf [Z ~> m]
  type(time_type) :: dummy_time
!
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  dummy_time = set_time(0,0)
  isd=G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    do i=isd,ied
      OD = CS%bed_elev(i,j) - rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf)
      if (OD >= 0) then
    ! ice thickness does not take up whole ocean column -> floating
        CS%OD_av(i,j) = OD
        CS%ground_frac(i,j) = 0.
      else
        CS%OD_av(i,j) = 0.
        CS%ground_frac(i,j) = 1.
      endif
    enddo
  enddo

  call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf,CS%taudx_shelf,CS%taudy_shelf, iters, Time)
end subroutine initialize_diagnostic_fields

!> This function returns the global maximum advective timestep that can be taken based on the current
!! ice velocities.  Because it involves finding a global minimum, it can be surprisingly expensive.
function ice_time_step_CFL(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G   !< The grid structure used by the ice shelf.
  real :: ice_time_step_CFL !< The maximum permitted timestep based on the ice velocities [T ~> s].

  real :: dt_local, min_dt ! These should be the minimum stable timesteps at a CFL of 1 [T ~> s]
  real :: min_vel          ! A minimal velocity for estimating a timestep [L T-1 ~> m s-1]
  integer :: i, j

  min_dt = 5.0e17*G%US%s_to_T ! The starting maximum is roughly the lifetime of the universe.
  min_vel = (1.0e-12/(365.0*86400.0)) * G%US%m_s_to_L_T
  do j=G%jsc,G%jec ; do i=G%isc,G%iec ; if (ISS%hmask(i,j) == 1.0 .or. ISS%hmask(i,j)==3) then
    dt_local = 2.0*G%areaT(i,j) / &
       (((G%dyCu(I,j)  * max(abs(CS%u_shelf(I,J)  + CS%u_shelf(I,j-1)), min_vel)) + &
         (G%dyCu(I-1,j)* max(abs(CS%u_shelf(I-1,J)+ CS%u_shelf(I-1,j-1)), min_vel))) + &
        ((G%dxCv(i,J)  * max(abs(CS%v_shelf(i,J)  + CS%v_shelf(i-1,J)), min_vel)) + &
         (G%dxCv(i,J-1)* max(abs(CS%v_shelf(i,J-1)+ CS%v_shelf(i-1,J-1)), min_vel))))

    min_dt = min(min_dt, dt_local)
  endif ; enddo ; enddo ! i- and j- loops

  call min_across_PEs(min_dt)

  ice_time_step_CFL = CS%CFL_factor * min_dt

end function ice_time_step_CFL

!> This subroutine updates the ice shelf velocities, mass, stresses and properties due to the
!! ice shelf dynamics.
subroutine update_ice_shelf(CS, ISS, G, US, time_step, Time, calve_ice_shelf_bergs, &
                            ocean_mass, coupled_grounding, must_update_vel)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                              !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time
  logical,                intent(in)    :: calve_ice_shelf_bergs !< To convert ice flux through front
                                                                 !! to bergs
  real, dimension(SZDI_(G),SZDJ_(G)), &
                optional, intent(in)    :: ocean_mass !< If present this is the mass per unit area
                                              !! of the ocean [R Z ~> kg m-2].
  logical,      optional, intent(in)    :: coupled_grounding !< If true, the grounding line is
                                              !! determined by coupled ice-ocean dynamics
  logical,      optional, intent(in)    :: must_update_vel !< Always update the ice velocities if true.
  integer :: iters
  logical :: update_ice_vel, coupled_GL

  update_ice_vel = .false.
  if (present(must_update_vel)) update_ice_vel = must_update_vel

  coupled_GL = .false.
  if (present(ocean_mass) .and. present(coupled_grounding)) coupled_GL = coupled_grounding
!
  if (CS%advect_shelf) then
    call ice_shelf_advect(CS, ISS, G, time_step, Time, calve_ice_shelf_bergs)
    if (CS%alternate_first_direction_IS) then
      CS%first_direction_IS = modulo(CS%first_direction_IS+1,2)
      CS%first_dir_restart_IS = real(CS%first_direction_IS)
    endif
  endif
  CS%elapsed_velocity_time = CS%elapsed_velocity_time + time_step
  if (CS%elapsed_velocity_time >= CS%velocity_update_time_step) update_ice_vel = .true.

  if (coupled_GL) then
    call update_OD_ffrac(CS, G, US, ocean_mass, update_ice_vel)
  elseif (update_ice_vel) then
    call update_OD_ffrac_uncoupled(CS, G, ISS%h_shelf(:,:))
    CS%GL_couple=.false.
  endif

  if (update_ice_vel) then
    call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf,CS%taudx_shelf,CS%taudy_shelf, iters, Time)
    CS%elapsed_velocity_time = 0.0
  endif

! call ice_shelf_temp(CS, ISS, G, US, time_step, ISS%water_flux, Time)

end subroutine update_ice_shelf

subroutine volume_above_floatation(CS, G, ISS, vaf, hemisphere)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  type(ice_shelf_state),  intent(in) :: ISS !< A structure with elements that describe
                                            !! the ice-shelf state
  real, intent(out) :: vaf !< area integrated volume above floatation [Z L2 ~> m3]
  integer, optional, intent(in) :: hemisphere !< 0 for Antarctica only, 1 for Greenland only. Otherwise, all ice sheets
  integer :: IS_ID ! local copy of hemisphere
  real, dimension(SZI_(G),SZJ_(G))  :: vaf_cell !< cell-wise volume above floatation [Z L2 ~> m3]
  integer, dimension(SZI_(G),SZJ_(G))  :: mask ! a mask for active cells depending on hemisphere indicated
  integer :: is,ie,js,je,i,j
  real :: rhoi_rhow, rhow_rhoi

  if (CS%GL_couple) &
    call MOM_error(FATAL, "MOM_ice_shelf_dyn, volume above floatation calculation assumes GL_couple=.FALSE..")

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  rhow_rhoi = CS%density_ocean_avg / CS%density_ice
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec

  if (present(hemisphere)) then
    IS_ID=hemisphere
  else
    IS_ID=-1
  endif

  mask(:,:)=0
  if (IS_ID==0) then     !Antarctica (S. Hemisphere) only
    do j = js,je ; do i = is,ie
      if (ISS%hmask(i,j)>0 .and. G%geoLatT(i,j)<=0.0) mask(i,j)=1
    enddo ; enddo
  elseif (IS_ID==1) then !Greenland (N. Hemisphere) only
    do j = js,je ; do i = is,ie
      if (ISS%hmask(i,j)>0 .and. G%geoLatT(i,j)>0.0)  mask(i,j)=1
    enddo ; enddo
  else                   !All ice sheets
    mask(is:ie,js:je)=ISS%hmask(is:ie,js:je)
  endif

  vaf_cell(:,:)=0.0
  do j = js,je ; do i = is,ie
    if (mask(i,j)>0) then
      if (CS%bed_elev(i,j) <= 0) then
        !grounded above sea level
        vaf_cell(i,j) = ISS%h_shelf(i,j) * ISS%area_shelf_h(i,j)
      else
        !grounded if vaf_cell(i,j) > 0
        vaf_cell(i,j) = max(ISS%h_shelf(i,j) - rhow_rhoi * CS%bed_elev(i,j), 0.0) * ISS%area_shelf_h(i,j)
      endif
    endif
  enddo ; enddo

  vaf = reproducing_sum(vaf_cell, unscale=G%US%Z_to_m*G%US%L_to_m**2)
end subroutine volume_above_floatation

!> multiplies a variable with the ice sheet grounding fraction
subroutine masked_var_grounded(G,CS,var,varout)
  type(ocean_grid_type), intent(in) :: G !< The grid structure used by the ice shelf.
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  real, dimension(SZI_(G),SZJ_(G)), intent(in)  :: var !< variable in
  real, dimension(SZI_(G),SZJ_(G)), intent(out)  :: varout !<variable out
  integer :: i, j
  do j = G%jsc,G%jec ; do i = G%isc,G%iec
      varout(i,j) = var(i,j) * CS%ground_frac(i,j)
  enddo ; enddo
end subroutine masked_var_grounded

!> Ice shelf dynamics post_data calls
subroutine IS_dynamics_post_data(time_step, Time, CS, ISS, G)
  real :: time_step !< Length of time for post data averaging [T ~> s].
  type(time_type),        intent(in)    :: Time !< The current model time
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJB_(G))  :: taud_x, taud_y, taud  ! area-averaged driving stress [R L2 T-2 ~> Pa]
  real, dimension(SZDI_(G),SZDJ_(G))  :: ice_visc ! area-averaged vertically integrated ice viscosity
                                                  !! [R L2 Z T-1 ~> Pa s m]
  real, dimension(SZDI_(G),SZDJ_(G))  :: basal_tr ! area-averaged taub_beta field related to basal traction,
                                                  !! [R L T-1 ~> Pa s m-1]
  real, dimension(SZDI_(G),SZDJ_(G))   :: surf_slope ! the surface slope of the ice shelf/sheet [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: ice_speed ! ice sheet flow speed [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_jump_n     ! max-min of co-located DG(1) corner thickness
                                                       !! across the up to 4 cells (hmask=1) touching each
                                                       !! B-grid node [Z ~> m], 0 where <2 touching cells
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_node_mx   ! per-B-node max of co-located DG(1) corner
                                                       ! thickness across touching hmask=1 cells [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_node_mn   ! per-B-node min of co-located DG(1) corner
                                                       ! thickness across touching hmask=1 cells [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_jump_n_rel ! h_jump_n normalised by the mean cell-mean
                                                       !! thickness over the same touching cells [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Hmax_Bd, Hmin_Bd ! per-B-node cell-mean envelope built from
                                                           !! hmask=1 and hmask=3 cells, identical to
                                                           !! the one used by the nodal AFC limiter
                                                           !! [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: count_Bd     ! number of contributing cells per B-node [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_jump_env   ! Hmax_B - Hmin_B at each B-node [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_jump_env_rel ! Envelope width / mean Hbar [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_over_n     ! BJ-overshoot magnitude at each B-node [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_over_n_rel ! Overshoot normalised by envelope mean [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Smax_Bd, Smin_Bd ! Per-B-node cell-mean surface envelope [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: s_over_n     ! BJ-overshoot magnitude in s at each B-node [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: s_over_n_rel ! s-overshoot normalised by envelope mean [nondim]
  real :: cell_mean_s_d                            ! Cell-mean surface contribution [Z ~> m]
  real :: s_cSW, s_cSE, s_cNW, s_cNE               ! Corner surface elevations at a B-node [Z ~> m]
  real :: bed_B                                    ! Bed elevation at the B-node, single-valued [Z ~> m]
  real :: over_b_s                                 ! Per-corner BJ violation in s at a B-node [Z ~> m]
  real :: over_b                                  ! Per-corner BJ-bound violation [Z ~> m]
  real :: h_cSW, h_cSE, h_cNW, h_cNE              ! corner thickness candidates at a B-node [Z ~> m]
  real :: Hb_SW, Hb_SE, Hb_NW, Hb_NE              ! cell-mean thickness for each touching cell [Z ~> m]
  real :: hmax_b, hmin_b                          ! max and min of valid corner candidates [Z ~> m]
  real :: Hbar_sum                                ! sum of valid cell-mean thicknesses [Z ~> m]
  real :: Hbar_avg                                ! arithmetic mean of valid cell-mean thicknesses [Z ~> m]
  real :: cell_mean_val_d                         ! envelope contribution from one cell [Z ~> m]
  real, parameter :: H_LARGE_D = 1.0e30           ! Sentinel for "no contributing cell"
  integer :: n_valid                              ! number of touching cells with hmask==1 [nondim]
  logical :: vSW, vSE, vNW, vNE                   ! per-touching-cell validity flags
  integer :: ii, jj                               ! touching-cell indices on the T-grid
  real, dimension(SZDIB_(G),SZDJ_(G)) :: hjump_fu ! per-u-face DG(1) thickness jump max|[h]| [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: sjump_fu ! per-u-face surface-elevation jump max|[s]| [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: sjump_fu_rel ! per-u-face |[s]|/Hbar_avg [nondim]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: hjump_fu_sgn ! signed u-face [h] = mean(h_plus - h_minus) [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: sjump_fu_sgn ! signed u-face [s] = mean(s_plus - s_minus) [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: un_fu    ! per-u-face |u.n| [L T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJB_(G)) :: hjump_fv ! per-v-face DG(1) thickness jump max|[h]| [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: sjump_fv ! per-v-face surface-elevation jump max|[s]| [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: sjump_fv_rel ! per-v-face |[s]|/Hbar_avg [nondim]
  real, dimension(SZDI_(G),SZDJB_(G)) :: hjump_fv_sgn ! signed v-face [h] = mean(h_plus - h_minus) [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: sjump_fv_sgn ! signed v-face [s] = mean(s_plus - s_minus) [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: un_fv    ! per-v-face |v.n| [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: eps_fu   ! per-u-face eps_e [T-1 ~> s-1]
  real, dimension(SZDI_(G),SZDJB_(G)) :: eps_fv   ! per-v-face eps_e [T-1 ~> s-1]
  real :: u_mn_d, v_mn_d, u_pl_d, v_pl_d ! Side-averaged velocities for eps_e [L T-1 ~> m s-1]
  real :: dudx_d, dudy_d, dvdx_d, dvdy_d ! Face-midpoint velocity gradients [T-1 ~> s-1]
  integer :: i_lo_d, i_hi_d, j_lo_d, j_hi_d ! Boundary-aware neighbour indices
  real :: Hbar_face_avg                           ! 0.5*(Hbar_A + Hbar_B) per face [Z ~> m]
  real :: rr                                      ! ice/ocean density ratio [nondim]
  real :: bed1, bed2                              ! bed elevation at the 2 face-endpoint nodes [Z ~> m]
  real :: h_m1, h_p1, h_m2, h_p2                  ! minus/plus side corner thickness at nodes 1,2 [Z ~> m]
  real :: s_m1, s_p1, s_m2, s_p2                  ! minus/plus side surface elevation at nodes 1,2 [Z ~> m]

  integer :: i, j

    call enable_averages(time_step, Time, CS%diag)
    if (CS%id_col_thick > 0) call post_data(CS%id_col_thick, CS%OD_av, CS%diag)
    if (CS%id_u_shelf > 0) call post_data(CS%id_u_shelf, CS%u_shelf, CS%diag)
    if (CS%id_v_shelf > 0) call post_data(CS%id_v_shelf, CS%v_shelf, CS%diag)
    if (CS%id_shelf_speed > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        ice_speed(I,J) = sqrt((CS%u_shelf(I,J)**2) + (CS%v_shelf(I,J)**2))
      enddo ; enddo
      call post_data(CS%id_shelf_speed, ice_speed, CS%diag)
    endif
!   if (CS%id_t_shelf > 0) call post_data(CS%id_t_shelf, CS%t_shelf, CS%diag)
    if (CS%id_taudx_shelf > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        taud_x(I,J) = CS%taudx_shelf(I,J)*G%IareaBu(I,J)
      enddo ; enddo
      call post_data(CS%id_taudx_shelf, taud_x, CS%diag)
    endif
    if (CS%id_taudy_shelf > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        taud_y(I,J) = CS%taudy_shelf(I,J)*G%IareaBu(I,J)
      enddo ; enddo
      call post_data(CS%id_taudy_shelf, taud_y, CS%diag)
    endif
    if (CS%id_taud_shelf > 0) then
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        taud(I,J) = sqrt((CS%taudx_shelf(I,J)**2)+(CS%taudy_shelf(I,J)**2))*G%IareaBu(I,J)
      enddo ; enddo
      call post_data(CS%id_taud_shelf, taud, CS%diag)
    endif
    if (CS%id_sx_shelf > 0) call post_data(CS%id_sx_shelf, CS%sx_shelf, CS%diag)
    if (CS%id_sy_shelf > 0) call post_data(CS%id_sy_shelf, CS%sy_shelf, CS%diag)
    if (CS%id_surf_slope_mag_shelf > 0) then
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        surf_slope(i,j) = sqrt((CS%sx_shelf(i,j)**2)+(CS%sy_shelf(i,j)**2))
      enddo ; enddo
      call post_data(CS%id_surf_slope_mag_shelf, surf_slope, CS%diag)
    endif
    if (CS%id_ground_frac > 0) call post_data(CS%id_ground_frac, CS%ground_frac, CS%diag)
    if (CS%id_OD_av >0) call post_data(CS%id_OD_av, CS%OD_av,CS%diag)
    if (CS%id_visc_shelf > 0) then
      call ice_visc_diag(CS,G,ice_visc)
      call post_data(CS%id_visc_shelf, ice_visc, CS%diag)
    endif
    if (CS%id_taub > 0) then
      call calc_shelf_taub(CS, ISS, G, basal_tr)
      call post_data(CS%id_taub, basal_tr, CS%diag)
    endif
    if (CS%id_bed_node > 0) call post_data(CS%id_bed_node, CS%bed_node, CS%diag)
    if (CS%id_h_nodal_SW > 0) call post_data(CS%id_h_nodal_SW, CS%h_nodal(:,:,1,1), CS%diag)
    if (CS%id_h_nodal_SE > 0) call post_data(CS%id_h_nodal_SE, CS%h_nodal(:,:,2,1), CS%diag)
    if (CS%id_h_nodal_NW > 0) call post_data(CS%id_h_nodal_NW, CS%h_nodal(:,:,1,2), CS%diag)
    if (CS%id_h_nodal_NE > 0) call post_data(CS%id_h_nodal_NE, CS%h_nodal(:,:,2,2), CS%diag)
    if ((CS%id_h_jump_node > 0 .or. CS%id_h_jump_node_rel > 0 .or. &
         CS%id_h_node_max > 0 .or. CS%id_h_node_min > 0) .and. associated(CS%h_nodal)) then
      call pass_corner_field(CS%h_nodal, G)
      h_jump_n(:,:)     = 0.0
      h_jump_n_rel(:,:) = 0.0
      h_node_mx(:,:)    = 0.0
      h_node_mn(:,:)    = 0.0
      ! At B-node (I,J) the up to 4 touching T-cells are, in fixed order,
      !   SW=(I,J), SE=(I+1,J), NW=(I,J+1), NE=(I+1,J+1).
      ! Each contributes its DG(1) corner that is co-located at (I,J):
      !   SW -> h_nodal(I,  J,  2,2)   SE -> h_nodal(I+1,J,  1,2)
      !   NW -> h_nodal(I,  J+1,2,1)   NE -> h_nodal(I+1,J+1,1,1)
      ! Fixed traversal order is required so the result is bitwise identical
      ! under any horizontal decomposition (halo cells produce the same value)
      ! and under a 90 deg rotation of the grid (the labelling rotates with i,j).
      do J=G%JscB,G%JecB ; do I=G%IscB,G%IecB
        ii = I   ; jj = J     ; vSW = (ISS%hmask(ii,jj) == 1.0)
        h_cSW = 0.0 ; Hb_SW = 0.0
        if (vSW) then ; h_cSW = CS%h_nodal(ii,jj,2,2) ; Hb_SW = ISS%h_shelf(ii,jj) ; endif
        ii = I+1 ; jj = J     ; vSE = (ISS%hmask(ii,jj) == 1.0)
        h_cSE = 0.0 ; Hb_SE = 0.0
        if (vSE) then ; h_cSE = CS%h_nodal(ii,jj,1,2) ; Hb_SE = ISS%h_shelf(ii,jj) ; endif
        ii = I   ; jj = J+1   ; vNW = (ISS%hmask(ii,jj) == 1.0)
        h_cNW = 0.0 ; Hb_NW = 0.0
        if (vNW) then ; h_cNW = CS%h_nodal(ii,jj,2,1) ; Hb_NW = ISS%h_shelf(ii,jj) ; endif
        ii = I+1 ; jj = J+1   ; vNE = (ISS%hmask(ii,jj) == 1.0)
        h_cNE = 0.0 ; Hb_NE = 0.0
        if (vNE) then ; h_cNE = CS%h_nodal(ii,jj,1,1) ; Hb_NE = ISS%h_shelf(ii,jj) ; endif
        n_valid = 0
        if (vSW) n_valid = n_valid + 1
        if (vSE) n_valid = n_valid + 1
        if (vNW) n_valid = n_valid + 1
        if (vNE) n_valid = n_valid + 1
        if (n_valid < 2) cycle
        ! Initialise min/max from the first valid candidate (SW->SE->NW->NE).
        if (vSW) then
          hmax_b = h_cSW ; hmin_b = h_cSW
        elseif (vSE) then
          hmax_b = h_cSE ; hmin_b = h_cSE
        elseif (vNW) then
          hmax_b = h_cNW ; hmin_b = h_cNW
        else
          hmax_b = h_cNE ; hmin_b = h_cNE
        endif
        if (vSE) then
          if (h_cSE > hmax_b) hmax_b = h_cSE
          if (h_cSE < hmin_b) hmin_b = h_cSE
        endif
        if (vNW) then
          if (h_cNW > hmax_b) hmax_b = h_cNW
          if (h_cNW < hmin_b) hmin_b = h_cNW
        endif
        if (vNE) then
          if (h_cNE > hmax_b) hmax_b = h_cNE
          if (h_cNE < hmin_b) hmin_b = h_cNE
        endif
        h_jump_n(I,J) = hmax_b - hmin_b
        h_node_mx(I,J) = hmax_b
        h_node_mn(I,J) = hmin_b
        ! Mean cell-mean thickness over the same valid cells, with a fixed
        ! summation order to avoid FMA / reassociation differences.
        Hbar_sum = (Hb_SW + Hb_SE) + (Hb_NW + Hb_NE)
        Hbar_avg = Hbar_sum / real(n_valid)
        h_jump_n_rel(I,J) = h_jump_n(I,J) / max(CS%min_h_shelf, Hbar_avg)
      enddo ; enddo
      if (CS%id_h_jump_node     > 0) call post_data(CS%id_h_jump_node,     h_jump_n,     CS%diag)
      if (CS%id_h_jump_node_rel > 0) call post_data(CS%id_h_jump_node_rel, h_jump_n_rel, CS%diag)
      if (CS%id_h_node_max      > 0) call post_data(CS%id_h_node_max,      h_node_mx,    CS%diag)
      if (CS%id_h_node_min      > 0) call post_data(CS%id_h_node_min,      h_node_mn,    CS%diag)
    endif
    if ((CS%id_h_jump_envelope > 0 .or. CS%id_h_jump_envelope_rel > 0 .or. &
         CS%id_h_overshoot_node > 0 .or. CS%id_h_overshoot_node_rel > 0 .or. &
         CS%id_s_overshoot_node > 0 .or. CS%id_s_overshoot_node_rel > 0) .and. &
        associated(CS%h_nodal)) then
      ! Per-B-node cell-mean envelope width Hmax_B - Hmin_B and (for s-overshoot)
      ! the analogous surface envelope [Smin_B, Smax_B], built from the cell-mean
      ! thickness projected with the cell-center bed under flotation.
      call pass_corner_field(CS%h_nodal, G)
      rr = CS%density_ice / CS%density_ocean_avg
      Hmax_Bd(:,:) = -H_LARGE_D
      Hmin_Bd(:,:) =  H_LARGE_D
      Smax_Bd(:,:) = -H_LARGE_D
      Smin_Bd(:,:) =  H_LARGE_D
      count_Bd(:,:) = 0.0
      do j = G%jsd, G%jed ; do i = G%isd, G%ied
        if (ISS%hmask(i,j) == 1.0) then
          cell_mean_val_d = nodal_cell_mean(CS%h_nodal(i,j,:,:), CS%cell_mean_w(i,j,:,:))
        elseif (ISS%hmask(i,j) == 3.0) then
          cell_mean_val_d = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
        else
          cycle
        endif
        if (rr*cell_mean_val_d - CS%bed_elev(i,j) > 0.0) then
          cell_mean_s_d = cell_mean_val_d - CS%bed_elev(i,j)
        else
          cell_mean_s_d = (1.0 - rr)*cell_mean_val_d
        endif
        if (i-1 >= G%IsdB .and. j-1 >= G%JsdB) then
          Hmax_Bd(i-1, j-1) = max(Hmax_Bd(i-1, j-1), cell_mean_val_d)
          Hmin_Bd(i-1, j-1) = min(Hmin_Bd(i-1, j-1), cell_mean_val_d)
          Smax_Bd(i-1, j-1) = max(Smax_Bd(i-1, j-1), cell_mean_s_d)
          Smin_Bd(i-1, j-1) = min(Smin_Bd(i-1, j-1), cell_mean_s_d)
          count_Bd(i-1, j-1) = count_Bd(i-1, j-1) + 1.0
        endif
        if (i <= G%IedB .and. j-1 >= G%JsdB) then
          Hmax_Bd(i,   j-1) = max(Hmax_Bd(i,   j-1), cell_mean_val_d)
          Hmin_Bd(i,   j-1) = min(Hmin_Bd(i,   j-1), cell_mean_val_d)
          Smax_Bd(i,   j-1) = max(Smax_Bd(i,   j-1), cell_mean_s_d)
          Smin_Bd(i,   j-1) = min(Smin_Bd(i,   j-1), cell_mean_s_d)
          count_Bd(i,   j-1) = count_Bd(i,   j-1) + 1.0
        endif
        if (i-1 >= G%IsdB .and. j <= G%JedB) then
          Hmax_Bd(i-1, j  ) = max(Hmax_Bd(i-1, j  ), cell_mean_val_d)
          Hmin_Bd(i-1, j  ) = min(Hmin_Bd(i-1, j  ), cell_mean_val_d)
          Smax_Bd(i-1, j  ) = max(Smax_Bd(i-1, j  ), cell_mean_s_d)
          Smin_Bd(i-1, j  ) = min(Smin_Bd(i-1, j  ), cell_mean_s_d)
          count_Bd(i-1, j  ) = count_Bd(i-1, j  ) + 1.0
        endif
        if (i <= G%IedB .and. j <= G%JedB) then
          Hmax_Bd(i,   j  ) = max(Hmax_Bd(i,   j  ), cell_mean_val_d)
          Hmin_Bd(i,   j  ) = min(Hmin_Bd(i,   j  ), cell_mean_val_d)
          Smax_Bd(i,   j  ) = max(Smax_Bd(i,   j  ), cell_mean_s_d)
          Smin_Bd(i,   j  ) = min(Smin_Bd(i,   j  ), cell_mean_s_d)
          count_Bd(i,   j  ) = count_Bd(i,   j  ) + 1.0
        endif
      enddo ; enddo
      call pass_var(Hmax_Bd,  G%domain, position=CORNER)
      call pass_var(Hmin_Bd,  G%domain, position=CORNER)
      call pass_var(Smax_Bd,  G%domain, position=CORNER)
      call pass_var(Smin_Bd,  G%domain, position=CORNER)
      call pass_var(count_Bd, G%domain, position=CORNER)
      h_jump_env(:,:)     = 0.0
      h_jump_env_rel(:,:) = 0.0
      do J = G%JscB, G%JecB ; do I = G%IscB, G%IecB
        if (count_Bd(I,J) < 1.5) cycle
        if (Hmax_Bd(I,J) <= -H_LARGE_D + 1.0 .or. Hmin_Bd(I,J) >= H_LARGE_D - 1.0) cycle
        h_jump_env(I,J)     = Hmax_Bd(I,J) - Hmin_Bd(I,J)
        h_jump_env_rel(I,J) = h_jump_env(I,J) / &
                              max(CS%min_h_shelf, 0.5*(Hmax_Bd(I,J) + Hmin_Bd(I,J)))
      enddo ; enddo
      if (CS%id_h_jump_envelope     > 0) call post_data(CS%id_h_jump_envelope,     h_jump_env,     CS%diag)
      if (CS%id_h_jump_envelope_rel > 0) call post_data(CS%id_h_jump_envelope_rel, h_jump_env_rel, CS%diag)
      ! BJ-overshoot per B-node: max over the up-to-4 touching DG(1) corner values of
      ! how far the corner sits outside the local cell-mean envelope [Hmin_B, Hmax_B].
      ! Nonzero only where the DG(1) corner has overshot the neighbor cell-mean
      ! envelope, i.e. the discontinuous mode is doing more than just resolve a sharp
      ! gradient between neighbouring cells. Uses Hmax_Bd / Hmin_Bd built above and
      ! the same corner-gathering pattern as h_jump_node.
      if (CS%id_h_overshoot_node > 0 .or. CS%id_h_overshoot_node_rel > 0) then
        h_over_n(:,:)     = 0.0
        h_over_n_rel(:,:) = 0.0
        do J = G%JscB, G%JecB ; do I = G%IscB, G%IecB
          if (count_Bd(I,J) < 1.5) cycle
          if (Hmax_Bd(I,J) <= -H_LARGE_D + 1.0 .or. Hmin_Bd(I,J) >= H_LARGE_D - 1.0) cycle
          ii = I   ; jj = J     ; vSW = (ISS%hmask(ii,jj) == 1.0)
          h_cSW = 0.0
          if (vSW) h_cSW = CS%h_nodal(ii,jj,2,2)
          ii = I+1 ; jj = J     ; vSE = (ISS%hmask(ii,jj) == 1.0)
          h_cSE = 0.0
          if (vSE) h_cSE = CS%h_nodal(ii,jj,1,2)
          ii = I   ; jj = J+1   ; vNW = (ISS%hmask(ii,jj) == 1.0)
          h_cNW = 0.0
          if (vNW) h_cNW = CS%h_nodal(ii,jj,2,1)
          ii = I+1 ; jj = J+1   ; vNE = (ISS%hmask(ii,jj) == 1.0)
          h_cNE = 0.0
          if (vNE) h_cNE = CS%h_nodal(ii,jj,1,1)
          over_b = 0.0
          if (vSW) over_b = max(over_b, h_cSW - Hmax_Bd(I,J), Hmin_Bd(I,J) - h_cSW)
          if (vSE) over_b = max(over_b, h_cSE - Hmax_Bd(I,J), Hmin_Bd(I,J) - h_cSE)
          if (vNW) over_b = max(over_b, h_cNW - Hmax_Bd(I,J), Hmin_Bd(I,J) - h_cNW)
          if (vNE) over_b = max(over_b, h_cNE - Hmax_Bd(I,J), Hmin_Bd(I,J) - h_cNE)
          h_over_n(I,J)     = over_b
          h_over_n_rel(I,J) = over_b / max(CS%min_h_shelf, 0.5*(Hmax_Bd(I,J) + Hmin_Bd(I,J)))
        enddo ; enddo
        if (CS%id_h_overshoot_node     > 0) call post_data(CS%id_h_overshoot_node,     h_over_n,     CS%diag)
        if (CS%id_h_overshoot_node_rel > 0) call post_data(CS%id_h_overshoot_node_rel, h_over_n_rel, CS%diag)
      endif
      ! BJ-overshoot in the surface-elevation field. Corner s computed by per-side
      ! flotation from h_corner and the B-node bed (single-valued at the node, so all
      ! co-located corners share bed_B), then compared against the cell-mean s envelope.
      ! Discriminator for the driving-stress-relevant spurious surface mode: small
      ! s-overshoot with large h-overshoot is the harmless flotation-only effect, large
      ! s-overshoot is a real spurious surface oscillation that corrupts driving stress.
      if (CS%id_s_overshoot_node > 0 .or. CS%id_s_overshoot_node_rel > 0) then
        call pass_var(CS%bed_node, G%domain, position=CORNER)
        s_over_n(:,:)     = 0.0
        s_over_n_rel(:,:) = 0.0
        do J = G%JscB, G%JecB ; do I = G%IscB, G%IecB
          if (count_Bd(I,J) < 1.5) cycle
          if (Smax_Bd(I,J) <= -H_LARGE_D + 1.0 .or. Smin_Bd(I,J) >= H_LARGE_D - 1.0) cycle
          bed_B = CS%bed_node(I,J)
          ii = I   ; jj = J     ; vSW = (ISS%hmask(ii,jj) == 1.0)
          s_cSW = 0.0
          if (vSW) then
            h_cSW = CS%h_nodal(ii,jj,2,2)
            s_cSW = merge(h_cSW - bed_B, (1.0 - rr)*h_cSW, rr*h_cSW - bed_B > 0.0)
          endif
          ii = I+1 ; jj = J     ; vSE = (ISS%hmask(ii,jj) == 1.0)
          s_cSE = 0.0
          if (vSE) then
            h_cSE = CS%h_nodal(ii,jj,1,2)
            s_cSE = merge(h_cSE - bed_B, (1.0 - rr)*h_cSE, rr*h_cSE - bed_B > 0.0)
          endif
          ii = I   ; jj = J+1   ; vNW = (ISS%hmask(ii,jj) == 1.0)
          s_cNW = 0.0
          if (vNW) then
            h_cNW = CS%h_nodal(ii,jj,2,1)
            s_cNW = merge(h_cNW - bed_B, (1.0 - rr)*h_cNW, rr*h_cNW - bed_B > 0.0)
          endif
          ii = I+1 ; jj = J+1   ; vNE = (ISS%hmask(ii,jj) == 1.0)
          s_cNE = 0.0
          if (vNE) then
            h_cNE = CS%h_nodal(ii,jj,1,1)
            s_cNE = merge(h_cNE - bed_B, (1.0 - rr)*h_cNE, rr*h_cNE - bed_B > 0.0)
          endif
          over_b_s = 0.0
          if (vSW) over_b_s = max(over_b_s, s_cSW - Smax_Bd(I,J), Smin_Bd(I,J) - s_cSW)
          if (vSE) over_b_s = max(over_b_s, s_cSE - Smax_Bd(I,J), Smin_Bd(I,J) - s_cSE)
          if (vNW) over_b_s = max(over_b_s, s_cNW - Smax_Bd(I,J), Smin_Bd(I,J) - s_cNW)
          if (vNE) over_b_s = max(over_b_s, s_cNE - Smax_Bd(I,J), Smin_Bd(I,J) - s_cNE)
          s_over_n(I,J)     = over_b_s
          s_over_n_rel(I,J) = over_b_s / max(CS%min_h_shelf, 0.5*(abs(Smax_Bd(I,J)) + abs(Smin_Bd(I,J))))
        enddo ; enddo
        if (CS%id_s_overshoot_node     > 0) call post_data(CS%id_s_overshoot_node,     s_over_n,     CS%diag)
        if (CS%id_s_overshoot_node_rel > 0) call post_data(CS%id_s_overshoot_node_rel, s_over_n_rel, CS%diag)
      endif
    endif
    if ((CS%id_h_jump_face_u > 0 .or. CS%id_h_jump_face_v > 0 .or. &
         CS%id_s_jump_face_u > 0 .or. CS%id_s_jump_face_v > 0 .or. &
         CS%id_s_jump_face_u_rel > 0 .or. CS%id_s_jump_face_v_rel > 0 .or. &
         CS%id_un_face_u > 0 .or. CS%id_un_face_v > 0) .and. associated(CS%h_nodal)) then
      ! Per-face inter-element jump of the broken-Q1 thickness, split onto u-faces
      ! (Cu) and v-faces (Cv) and reported in both thickness [h] and surface
      ! elevation [s]. Unlike h_jump_node (max-min over the up-to-4 corners at a
      ! B-node), this attributes the discontinuity to a single face, so it can be
      ! correlated against un_face to test whether jumps accumulate at shear-margin
      ! faces where u.n ~ 0 (the advectively-uncoupled, undamped jump mode).
      call pass_corner_field(CS%h_nodal, G)
      call pass_var(CS%bed_node, G%domain, position=CORNER)
      call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE)
      rr = CS%density_ice / CS%density_ocean_avg
      hjump_fu(:,:) = 0.0 ; sjump_fu(:,:) = 0.0 ; sjump_fu_rel(:,:) = 0.0 ; un_fu(:,:) = 0.0
      hjump_fu_sgn(:,:) = 0.0 ; sjump_fu_sgn(:,:) = 0.0
      hjump_fv(:,:) = 0.0 ; sjump_fv(:,:) = 0.0 ; sjump_fv_rel(:,:) = 0.0 ; un_fv(:,:) = 0.0
      hjump_fv_sgn(:,:) = 0.0 ; sjump_fv_sgn(:,:) = 0.0
      ! u-faces: minus side = west cell (I,j) east edge, plus side = east cell
      ! (I+1,j) west edge; endpoint nodes 1=south (I,j-1), 2=north (I,j).
      do j = G%jsc, G%jec ; do I = G%IscB, G%IecB
        if (ISS%hmask(I,j) /= 1.0) cycle
        if (ISS%hmask(I+1,j) /= 1.0) cycle
        bed1 = CS%bed_node(I,j-1) ; bed2 = CS%bed_node(I,j)
        h_m1 = CS%h_nodal(I,  j,2,1) ; h_p1 = CS%h_nodal(I+1,j,1,1)
        h_m2 = CS%h_nodal(I,  j,2,2) ; h_p2 = CS%h_nodal(I+1,j,1,2)
        s_m1 = merge(h_m1-bed1, (1.0-rr)*h_m1, rr*h_m1-bed1 > 0.0)
        s_p1 = merge(h_p1-bed1, (1.0-rr)*h_p1, rr*h_p1-bed1 > 0.0)
        s_m2 = merge(h_m2-bed2, (1.0-rr)*h_m2, rr*h_m2-bed2 > 0.0)
        s_p2 = merge(h_p2-bed2, (1.0-rr)*h_p2, rr*h_p2-bed2 > 0.0)
        hjump_fu(I,j) = max(abs(h_m1-h_p1), abs(h_m2-h_p2))
        sjump_fu(I,j) = max(abs(s_m1-s_p1), abs(s_m2-s_p2))
        hjump_fu_sgn(I,j) = 0.5*((h_p1 - h_m1) + (h_p2 - h_m2))
        sjump_fu_sgn(I,j) = 0.5*((s_p1 - s_m1) + (s_p2 - s_m2))
        Hbar_face_avg = max(CS%min_h_shelf, 0.125*( &
          (CS%h_nodal(I  ,j,1,1) + CS%h_nodal(I  ,j,2,1)) + &
          (CS%h_nodal(I  ,j,1,2) + CS%h_nodal(I  ,j,2,2)) + &
          (CS%h_nodal(I+1,j,1,1) + CS%h_nodal(I+1,j,2,1)) + &
          (CS%h_nodal(I+1,j,1,2) + CS%h_nodal(I+1,j,2,2))))
        sjump_fu_rel(I,j) = sjump_fu(I,j) / Hbar_face_avg
        un_fu(I,j) = abs(0.5*(CS%u_shelf(I,j-1) + CS%u_shelf(I,j)))
      enddo ; enddo
      ! v-faces: minus side = south cell (i,J) north edge, plus side = north cell
      ! (i,J+1) south edge; endpoint nodes 1=west (i-1,J), 2=east (i,J).
      do J = G%JscB, G%JecB ; do i = G%isc, G%iec
        if (ISS%hmask(i,J) /= 1.0) cycle
        if (ISS%hmask(i,J+1) /= 1.0) cycle
        bed1 = CS%bed_node(i-1,J) ; bed2 = CS%bed_node(i,J)
        h_m1 = CS%h_nodal(i,J,  1,2) ; h_p1 = CS%h_nodal(i,J+1,1,1)
        h_m2 = CS%h_nodal(i,J,  2,2) ; h_p2 = CS%h_nodal(i,J+1,2,1)
        s_m1 = merge(h_m1-bed1, (1.0-rr)*h_m1, rr*h_m1-bed1 > 0.0)
        s_p1 = merge(h_p1-bed1, (1.0-rr)*h_p1, rr*h_p1-bed1 > 0.0)
        s_m2 = merge(h_m2-bed2, (1.0-rr)*h_m2, rr*h_m2-bed2 > 0.0)
        s_p2 = merge(h_p2-bed2, (1.0-rr)*h_p2, rr*h_p2-bed2 > 0.0)
        hjump_fv(i,J) = max(abs(h_m1-h_p1), abs(h_m2-h_p2))
        sjump_fv(i,J) = max(abs(s_m1-s_p1), abs(s_m2-s_p2))
        hjump_fv_sgn(i,J) = 0.5*((h_p1 - h_m1) + (h_p2 - h_m2))
        sjump_fv_sgn(i,J) = 0.5*((s_p1 - s_m1) + (s_p2 - s_m2))
        Hbar_face_avg = max(CS%min_h_shelf, 0.125*( &
          (CS%h_nodal(i,J  ,1,1) + CS%h_nodal(i,J  ,2,1)) + &
          (CS%h_nodal(i,J  ,1,2) + CS%h_nodal(i,J  ,2,2)) + &
          (CS%h_nodal(i,J+1,1,1) + CS%h_nodal(i,J+1,2,1)) + &
          (CS%h_nodal(i,J+1,1,2) + CS%h_nodal(i,J+1,2,2))))
        sjump_fv_rel(i,J) = sjump_fv(i,J) / Hbar_face_avg
        un_fv(i,J) = abs(0.5*(CS%v_shelf(i-1,J) + CS%v_shelf(i,J)))
      enddo ; enddo
      if (CS%id_h_jump_face_u > 0) call post_data(CS%id_h_jump_face_u, hjump_fu, CS%diag)
      if (CS%id_h_jump_face_v > 0) call post_data(CS%id_h_jump_face_v, hjump_fv, CS%diag)
      if (CS%id_s_jump_face_u > 0) call post_data(CS%id_s_jump_face_u, sjump_fu, CS%diag)
      if (CS%id_s_jump_face_v > 0) call post_data(CS%id_s_jump_face_v, sjump_fv, CS%diag)
      if (CS%id_h_jump_face_u_signed > 0) &
          call post_data(CS%id_h_jump_face_u_signed, hjump_fu_sgn, CS%diag)
      if (CS%id_h_jump_face_v_signed > 0) &
          call post_data(CS%id_h_jump_face_v_signed, hjump_fv_sgn, CS%diag)
      if (CS%id_s_jump_face_u_signed > 0) &
          call post_data(CS%id_s_jump_face_u_signed, sjump_fu_sgn, CS%diag)
      if (CS%id_s_jump_face_v_signed > 0) &
          call post_data(CS%id_s_jump_face_v_signed, sjump_fv_sgn, CS%diag)
      if (CS%id_s_jump_face_u_rel > 0) &
          call post_data(CS%id_s_jump_face_u_rel, sjump_fu_rel, CS%diag)
      if (CS%id_s_jump_face_v_rel > 0) &
          call post_data(CS%id_s_jump_face_v_rel, sjump_fv_rel, CS%diag)
      if (CS%id_un_face_u > 0) call post_data(CS%id_un_face_u, un_fu, CS%diag)
      if (CS%id_un_face_v > 0) call post_data(CS%id_un_face_v, un_fv, CS%diag)
    endif
    ! Per-face SSA effective strain rate eps_e using the same boundary-aware stencil
    ! as the DG(1) artificial viscosity. Posted unconditionally of art_visc on/off so
    ! it can be used to tune DG1_ART_VISC_STRAIN_COEF in advance.
    if (CS%id_dg_eps_face_u > 0) then
      eps_fu(:,:) = 0.0
      do j = G%jsc, G%jec ; do I = G%IscB, G%IecB
        if (ISS%hmask(I,  j) /= 1.0 .and. ISS%hmask(I,  j) /= 3.0) cycle
        if (ISS%hmask(I+1,j) /= 1.0 .and. ISS%hmask(I+1,j) /= 3.0) cycle
        i_lo_d = max(I-1, G%isd) ; i_hi_d = min(I+1, G%ied)
        if (i_lo_d < I) then
          u_mn_d = 0.25*((CS%u_shelf(i_lo_d,j-1) + CS%u_shelf(I,j-1)) + &
                         (CS%u_shelf(i_lo_d,j  ) + CS%u_shelf(I,j  )))
          v_mn_d = 0.25*((CS%v_shelf(i_lo_d,j-1) + CS%v_shelf(I,j-1)) + &
                         (CS%v_shelf(i_lo_d,j  ) + CS%v_shelf(I,j  )))
        else
          u_mn_d = 0.5*(CS%u_shelf(I,j-1) + CS%u_shelf(I,j))
          v_mn_d = 0.5*(CS%v_shelf(I,j-1) + CS%v_shelf(I,j))
        endif
        if (i_hi_d > I) then
          u_pl_d = 0.25*((CS%u_shelf(I    ,j-1) + CS%u_shelf(i_hi_d,j-1)) + &
                         (CS%u_shelf(I    ,j  ) + CS%u_shelf(i_hi_d,j  )))
          v_pl_d = 0.25*((CS%v_shelf(I    ,j-1) + CS%v_shelf(i_hi_d,j-1)) + &
                         (CS%v_shelf(I    ,j  ) + CS%v_shelf(i_hi_d,j  )))
        else
          u_pl_d = 0.5*(CS%u_shelf(I,j-1) + CS%u_shelf(I,j))
          v_pl_d = 0.5*(CS%v_shelf(I,j-1) + CS%v_shelf(I,j))
        endif
        dudx_d = (u_pl_d - u_mn_d) / G%dxCu(I,j)
        dvdx_d = (v_pl_d - v_mn_d) / G%dxCu(I,j)
        dudy_d = (CS%u_shelf(I,j) - CS%u_shelf(I,j-1)) / G%dyCu(I,j)
        dvdy_d = (CS%v_shelf(I,j) - CS%v_shelf(I,j-1)) / G%dyCu(I,j)
        eps_fu(I,j) = dg1_face_eps_eff(dudx_d, dudy_d, dvdx_d, dvdy_d)
      enddo ; enddo
      call post_data(CS%id_dg_eps_face_u, eps_fu, CS%diag)
    endif
    if (CS%id_dg_eps_face_v > 0) then
      eps_fv(:,:) = 0.0
      do J = G%JscB, G%JecB ; do i = G%isc, G%iec
        if (ISS%hmask(i,J  ) /= 1.0 .and. ISS%hmask(i,J  ) /= 3.0) cycle
        if (ISS%hmask(i,J+1) /= 1.0 .and. ISS%hmask(i,J+1) /= 3.0) cycle
        j_lo_d = max(J-1, G%jsd) ; j_hi_d = min(J+1, G%jed)
        if (j_lo_d < J) then
          u_mn_d = 0.25*((CS%u_shelf(i-1,j_lo_d) + CS%u_shelf(i,j_lo_d)) + &
                         (CS%u_shelf(i-1,J     ) + CS%u_shelf(i,J     )))
          v_mn_d = 0.25*((CS%v_shelf(i-1,j_lo_d) + CS%v_shelf(i,j_lo_d)) + &
                         (CS%v_shelf(i-1,J     ) + CS%v_shelf(i,J     )))
        else
          u_mn_d = 0.5*(CS%u_shelf(i-1,J) + CS%u_shelf(i,J))
          v_mn_d = 0.5*(CS%v_shelf(i-1,J) + CS%v_shelf(i,J))
        endif
        if (j_hi_d > J) then
          u_pl_d = 0.25*((CS%u_shelf(i-1,J     ) + CS%u_shelf(i,J     )) + &
                         (CS%u_shelf(i-1,j_hi_d) + CS%u_shelf(i,j_hi_d)))
          v_pl_d = 0.25*((CS%v_shelf(i-1,J     ) + CS%v_shelf(i,J     )) + &
                         (CS%v_shelf(i-1,j_hi_d) + CS%v_shelf(i,j_hi_d)))
        else
          u_pl_d = 0.5*(CS%u_shelf(i-1,J) + CS%u_shelf(i,J))
          v_pl_d = 0.5*(CS%v_shelf(i-1,J) + CS%v_shelf(i,J))
        endif
        dudy_d = (u_pl_d - u_mn_d) / G%dyCv(i,J)
        dvdy_d = (v_pl_d - v_mn_d) / G%dyCv(i,J)
        dudx_d = (CS%u_shelf(i,J) - CS%u_shelf(i-1,J)) / G%dxCv(i,J)
        dvdx_d = (CS%v_shelf(i,J) - CS%v_shelf(i-1,J)) / G%dxCv(i,J)
        eps_fv(i,J) = dg1_face_eps_eff(dudx_d, dudy_d, dvdx_d, dvdy_d)
      enddo ; enddo
      call post_data(CS%id_dg_eps_face_v, eps_fv, CS%diag)
    endif
    if (CS%id_h_source_rate > 0 .and. associated(CS%h_source_rate_last)) &
        call post_data(CS%id_h_source_rate, CS%h_source_rate_last, CS%diag)
    if (CS%id_phi_x_FV > 0 .and. associated(CS%phi_x_FV)) &
        call post_data(CS%id_phi_x_FV, CS%phi_x_FV, CS%diag)
    if (CS%id_phi_y_FV > 0 .and. associated(CS%phi_y_FV)) &
        call post_data(CS%id_phi_y_FV, CS%phi_y_FV, CS%diag)
    if (CS%id_dg_lim_phi_xi > 0 .and. associated(CS%dg_lim_phi_xi)) &
        call post_data(CS%id_dg_lim_phi_xi, CS%dg_lim_phi_xi, CS%diag)
    if (CS%id_dg_lim_phi_eta > 0 .and. associated(CS%dg_lim_phi_eta)) &
        call post_data(CS%id_dg_lim_phi_eta, CS%dg_lim_phi_eta, CS%diag)
    if (CS%id_dg_lim_phi_cross > 0 .and. associated(CS%dg_lim_phi_cross)) &
        call post_data(CS%id_dg_lim_phi_cross, CS%dg_lim_phi_cross, CS%diag)
    if (CS%id_dg_lim_mass_drift > 0 .and. associated(CS%dg_lim_mass_drift)) &
        call post_data(CS%id_dg_lim_mass_drift, CS%dg_lim_mass_drift, CS%diag)
    if (CS%id_dg_lim_phi > 0 .and. associated(CS%dg_lim_phi)) &
        call post_data(CS%id_dg_lim_phi, CS%dg_lim_phi, CS%diag)
    if (CS%id_dg_lim_pk_factor > 0 .and. associated(CS%dg_lim_pk_factor)) &
        call post_data(CS%id_dg_lim_pk_factor, CS%dg_lim_pk_factor, CS%diag)
    if (CS%id_dg_art_visc_coef_u > 0 .and. associated(CS%dg_art_visc_coef_u)) &
        call post_data(CS%id_dg_art_visc_coef_u, CS%dg_art_visc_coef_u, CS%diag)
    if (CS%id_dg_art_visc_coef_v > 0 .and. associated(CS%dg_art_visc_coef_v)) &
        call post_data(CS%id_dg_art_visc_coef_v, CS%dg_art_visc_coef_v, CS%diag)
    if (CS%id_dg_art_visc_nu_u > 0 .and. associated(CS%dg_art_visc_nu_u)) &
        call post_data(CS%id_dg_art_visc_nu_u, CS%dg_art_visc_nu_u, CS%diag)
    if (CS%id_dg_art_visc_nu_v > 0 .and. associated(CS%dg_art_visc_nu_v)) &
        call post_data(CS%id_dg_art_visc_nu_v, CS%dg_art_visc_nu_v, CS%diag)
    if (CS%id_dg_art_visc_excess_frac_u > 0 .and. associated(CS%dg_art_visc_excess_frac_u)) &
        call post_data(CS%id_dg_art_visc_excess_frac_u, CS%dg_art_visc_excess_frac_u, CS%diag)
    if (CS%id_dg_art_visc_excess_frac_v > 0 .and. associated(CS%dg_art_visc_excess_frac_v)) &
        call post_data(CS%id_dg_art_visc_excess_frac_v, CS%dg_art_visc_excess_frac_v, CS%diag)
    if (CS%id_dg_art_visc_allow_u > 0 .and. associated(CS%dg_art_visc_allow_u)) &
        call post_data(CS%id_dg_art_visc_allow_u, CS%dg_art_visc_allow_u, CS%diag)
    if (CS%id_dg_art_visc_allow_v > 0 .and. associated(CS%dg_art_visc_allow_v)) &
        call post_data(CS%id_dg_art_visc_allow_v, CS%dg_art_visc_allow_v, CS%diag)
    if (CS%id_dg_art_visc_cell_scale > 0 .and. associated(CS%dg_art_visc_cell_scale)) &
        call post_data(CS%id_dg_art_visc_cell_scale, CS%dg_art_visc_cell_scale, CS%diag)
    if (CS%id_dg_slow_idle_face_u > 0 .and. associated(CS%dg_slow_idle_face_u)) &
        call post_data(CS%id_dg_slow_idle_face_u, CS%dg_slow_idle_face_u, CS%diag)
    if (CS%id_dg_slow_idle_face_v > 0 .and. associated(CS%dg_slow_idle_face_v)) &
        call post_data(CS%id_dg_slow_idle_face_v, CS%dg_slow_idle_face_v, CS%diag)
    if (CS%id_u_mask > 0) call post_data(CS%id_u_mask, CS%umask, CS%diag)
    if (CS%id_v_mask > 0) call post_data(CS%id_v_mask, CS%vmask, CS%diag)
    if (CS%id_ufb_mask > 0) call post_data(CS%id_ufb_mask, CS%u_face_mask_bdry, CS%diag)
    if (CS%id_vfb_mask > 0) call post_data(CS%id_vfb_mask, CS%v_face_mask_bdry, CS%diag)
!   if (CS%id_t_mask > 0) call post_data(CS%id_t_mask, CS%tmask, CS%diag)

    if (CS%id_duHdx > 0         .or. CS%id_dvHdy > 0         .or. CS%id_fluxdiv > 0       .or. &
        CS%id_devstress_xx > 0  .or. CS%id_devstress_yy > 0  .or. CS%id_devstress_xy > 0  .or. &
        CS%id_strainrate_xx > 0 .or. CS%id_strainrate_yy > 0 .or. CS%id_strainrate_xy > 0 .or. &
        CS%id_pdevstress_1 > 0  .or. CS%id_pdevstress_2 > 0  .or. &
        CS%id_pstrainrate_1 > 0 .or. CS%id_pstrainrate_2 > 0) then
      call IS_dynamics_post_data_2(CS, ISS, G)
    endif

    call disable_averaging(CS%diag)
end subroutine IS_dynamics_post_data

!> Calculate cell-centered, area-averaged, vertically integrated ice viscosity for diagnostics
subroutine ice_visc_diag(CS,G,ice_visc)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(out)  :: ice_visc !< area-averaged vertically integrated ice viscosity
                                                               !! [R L2 Z T-1 ~> Pa s m]
  integer :: i,j

  ice_visc(:,:)=0.0
  if (CS%visc_qps==4) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ice_visc(i,j) = (0.25 * G%IareaT(i,j)) * &
        ((CS%ice_visc(i,j,1) + CS%ice_visc(i,j,4)) + (CS%ice_visc(i,j,2) + CS%ice_visc(i,j,3)))
    enddo ; enddo
  else
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      ice_visc(i,j) = CS%ice_visc(i,j,1)*G%IareaT(i,j)
    enddo ; enddo
  endif
end subroutine ice_visc_diag

!>  Writes the total ice shelf kinetic energy and mass to an ascii file
subroutine write_ice_shelf_energy(CS, G, US, mass, area, day, time_step, mass_hole)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: mass !< The mass per unit area of the ice shelf
                                                !! or sheet [R Z ~> kg m-2]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                           intent(in)    :: area !< The ice shelf or ice sheet area [L2 ~> m2]
  type(time_type),         intent(in)    :: day !< The current model time.
  type(time_type),  optional, intent(in) :: time_step !< The current time step
  real, optional, intent(in) :: mass_hole !< ice-sheet mass in the ocean grid hole, if present [RZL2 ~> kg]
  ! Local variables
  type(time_type) :: dt ! A time_type version of the timestep.
  real, dimension(SZDI_(G),SZDJ_(G)) :: tmp1 ! A temporary array used in reproducing sums [various]
  real :: KE_tot    ! The total kinetic energy [R Z L4 T-2 ~> J]
  real :: mass_tot  ! The total mass [R Z L2 ~> kg]
  integer :: is, ie, js, je, isr, ier, jsr, jer, i, j
  character(len=32)  :: mesg_intro, time_units, day_str, n_str, date_str
  integer :: start_of_day, num_days
  real    :: reday  ! Time in units given by CS%Timeunit, but often [days]

  ! write_energy_time is the next integral multiple of energysavedays.
  if (present(time_step)) then
    dt = time_step
  else
    dt = set_time(seconds=2)
  endif

   !CS%prev_IS_energy_calls tracks the ice sheet step, which is outputted in the energy file.
  if (CS%prev_IS_energy_calls == 0) then
    if (CS%energysave_geometric) then
      if (CS%energysavedays_geometric < CS%energysavedays) then
        CS%write_energy_time = day + CS%energysavedays_geometric
        CS%geometric_end_time = CS%Start_time + CS%energysavedays * &
          (1 + (day - CS%Start_time) / CS%energysavedays)
      else
        CS%write_energy_time = CS%Start_time + CS%energysavedays * &
          (1 + (day - CS%Start_time) / CS%energysavedays)
      endif
    else
      CS%write_energy_time = CS%Start_time + CS%energysavedays * &
        (1 + (day - CS%Start_time) / CS%energysavedays)
    endif
  elseif (day + (dt/2) <= CS%write_energy_time) then
    CS%prev_IS_energy_calls = CS%prev_IS_energy_calls + 1
    return  ! Do not write this step
  else ! Determine the next write time before proceeding
    if (CS%energysave_geometric) then
      if (CS%write_energy_time + CS%energysavedays_geometric >= &
          CS%geometric_end_time) then
        CS%write_energy_time = CS%geometric_end_time
        CS%energysave_geometric = .false.  ! stop geometric progression
      else
        CS%write_energy_time = CS%write_energy_time + CS%energysavedays_geometric
      endif
      CS%energysavedays_geometric = CS%energysavedays_geometric*2
    else
      CS%write_energy_time = CS%write_energy_time + CS%energysavedays
    endif
  endif

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isr = is - (G%isd-1) ; ier = ie - (G%isd-1) ; jsr = js - (G%jsd-1) ; jer = je - (G%jsd-1)

  !calculate KE using cell-centered ice shelf velocity
  tmp1(:,:) = 0.0
  do j=js,je ; do i=is,ie
    tmp1(i,j) = 0.03125 * (mass(i,j) * area(i,j)) * &
      ((((CS%u_shelf(I-1,J-1)+CS%u_shelf(I,J))+(CS%u_shelf(I,J-1)+CS%u_shelf(I-1,J)))**2) + &
       (((CS%v_shelf(I-1,J-1)+CS%v_shelf(I,J))+(CS%v_shelf(I,J-1)+CS%v_shelf(I-1,J)))**2))
  enddo ; enddo

  KE_tot = reproducing_sum(tmp1, isr, ier, jsr, jer, unscale=(US%RZL2_to_kg*US%L_T_to_m_s**2))

  !calculate mass
  tmp1(:,:) = 0.0
  do j=js,je ; do i=is,ie
    tmp1(i,j) = mass(i,j) * area(i,j)
  enddo ; enddo

  mass_tot = reproducing_sum(tmp1, isr, ier, jsr, jer, unscale=US%RZL2_to_kg)
  if (present(mass_hole)) mass_tot = mass_tot + mass_hole

  if (is_root_pe()) then  ! Only the root PE actually writes anything.
    if (day > CS%Start_time) then
      call open_ASCII_file(CS%IS_fileenergy_ascii, trim(CS%IS_energyfile), action=APPEND_FILE)
    else
      call open_ASCII_file(CS%IS_fileenergy_ascii, trim(CS%IS_energyfile), action=WRITEONLY_FILE)
      if (abs(CS%timeunit - 86400.0) < 1.0) then
        write(CS%IS_fileenergy_ascii,'("  Step,",7x,"Day,",8x,"Energy/Mass,",13x,"Total Mass")')
        write(CS%IS_fileenergy_ascii,'(12x,"[days]",10x,"[m2 s-2]",17x,"[kg]")')
      else
        if ((CS%timeunit >= 0.99) .and. (CS%timeunit < 1.01)) then
          time_units = "           [seconds]     "
        elseif ((CS%timeunit >= 3599.0) .and. (CS%timeunit < 3601.0)) then
          time_units = "            [hours]      "
        elseif ((CS%timeunit >= 86399.0) .and. (CS%timeunit < 86401.0)) then
          time_units = "             [days]      "
        elseif ((CS%timeunit >= 3.0e7) .and. (CS%timeunit < 3.2e7)) then
          time_units = "            [years]      "
        else
          write(time_units,'(9x,"[",es8.2," s]    ")') CS%timeunit
        endif

        write(CS%IS_fileenergy_ascii,'("  Step,",7x,"Time,",7x,"Energy/Mass,",13x,"Total Mass")')
        write(CS%IS_fileenergy_ascii,'(A25,3x,"[m2 s-2]",17x,"[kg]")') time_units
      endif
    endif

    call get_time(day, start_of_day, num_days)

    if (abs(CS%timeunit - 86400.0) < 1.0) then
      reday = REAL(num_days)+ (REAL(start_of_day)/86400.0)
    else
      reday = REAL(num_days)*(86400.0/CS%timeunit) + REAL(start_of_day)/abs(CS%timeunit)
    endif

    if (reday < 1.0e8) then ;      write(day_str, '(F12.3)') reday
    elseif (reday < 1.0e11) then ; write(day_str, '(F15.3)') reday
    else ;                         write(day_str, '(ES15.9)') reday ; endif

    if (CS%prev_IS_energy_calls < 1000000) then ; write(n_str, '(I6)') CS%prev_IS_energy_calls
    else ; write(n_str, '(I0)') CS%prev_IS_energy_calls ; endif

    write(CS%IS_fileenergy_ascii,'(A,",",A,", En ",ES22.16,", M ",ES11.5)') &
      trim(n_str), trim(day_str), US%L_T_to_m_s**2*KE_tot/mass_tot, US%RZL2_to_kg*mass_tot
  endif

  CS%prev_IS_energy_calls = CS%prev_IS_energy_calls + 1
end subroutine write_ice_shelf_energy

!> This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!! Additionally, it will update the volume of ice in partially-filled cells, and update
!! hmask accordingly
subroutine ice_shelf_advect(CS, ISS, G, time_step, Time, calve_ice_shelf_bergs)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< time step [T ~> s]
  type(time_type),        intent(in)    :: Time !< The current model time
  logical,                intent(in)    :: calve_ice_shelf_bergs !< If true, track ice shelf flux through a
                                               !! static ice shelf, so that it can be converted into icebergs

! 3/8/11 DNG
!
!    This subroutine takes the velocity (on the Bgrid) and timesteps h_t = - div (uh) once.
!    ADDITIONALLY, it will update the volume of ice in partially-filled cells, and update
!        hmask accordingly
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  real, dimension(SZDI_(G),SZDJ_(G))   :: h_after_flux1, h_after_flux2 ! Ice thicknesses [Z ~> m].
  real, dimension(SZDIB_(G),SZDJ_(G))  :: uh_ice  ! The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G))  :: vh_ice  ! The accumulated meridional ice volume flux [Z L2 ~> m3]
  type(loop_bounds_type) :: LB
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec, stencil

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  uh_ice(:,:) = 0.0
  vh_ice(:,:) = 0.0

  h_after_flux1(:,:) = 0.0
  h_after_flux2(:,:) = 0.0
  ! call MOM_mesg("MOM_ice_shelf.F90: ice_shelf_advect called")

  do j=jsd,jed ; do i=isd,ied ; if (CS%h_bdry_val(i,j) /= 0.0) then
    ISS%h_shelf(i,j) = CS%h_bdry_val(i,j)
  endif ; enddo ; enddo

  if (CS%use_DG_thickness) then
    ! Nodal DG(1) unsplit advection with SSP-RK2 on CS%h_nodal.
    call ice_shelf_advect_DG1_nodal(CS, ISS, G, time_step, ISS%hmask, uh_ice, vh_ice)
    ! Publish derived state: ISS%h_shelf as the cell-mean of CS%h_nodal.
    call recompute_h_shelf_from_nodal(CS, ISS, G)
    call pass_var(ISS%h_shelf, G%domain)
  else
    ! Reset FV slope-limiter face diagnostic; advect_thickness_x/y writes slope_lim
    ! at each face where the limiter branch is taken. Faces not visited (incomplete
    ! stencil at front, inactive face) report 1.0.
    if (associated(CS%phi_x_FV)) CS%phi_x_FV(:,:) = 1.0
    if (associated(CS%phi_y_FV)) CS%phi_y_FV(:,:) = 1.0
    stencil = 2
    if (modulo(CS%first_direction_IS,2)==0) then
      !x first
      LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc-stencil ; LB%jeh = G%jec+stencil
      if (LB%jsh < jsd) call MOM_error(FATAL, &
        "ice_shelf_advect:  Halo is too small for the ice thickness advection stencil.")
      call ice_shelf_advect_thickness_x(CS, G, LB, time_step, ISS%hmask, ISS%h_shelf, h_after_flux1, uh_ice)
      call pass_var(h_after_flux1, G%domain)
      LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
      call ice_shelf_advect_thickness_y(CS, G, LB, time_step, ISS%hmask, h_after_flux1, h_after_flux2, vh_ice)
    else
      ! y first
      LB%ish = G%isc-stencil ; LB%ieh = G%iec+stencil ; LB%jsh = G%jsc ; LB%jeh = G%jec
      if (LB%ish < isd) call MOM_error(FATAL, &
        "ice_shelf_advect:  Halo is too small for the ice thickness advection stencil.")
      call ice_shelf_advect_thickness_y(CS, G, LB, time_step, ISS%hmask, ISS%h_shelf, h_after_flux1, vh_ice)
      call pass_var(h_after_flux1, G%domain)
      LB%ish = G%isc ; LB%ieh = G%iec ; LB%jsh = G%jsc ; LB%jeh = G%jec
      call ice_shelf_advect_thickness_x(CS, G, LB, time_step, ISS%hmask, h_after_flux1, h_after_flux2, uh_ice)
    endif
    call pass_var(h_after_flux2, G%domain)

    do j=jsd,jed
      do i=isd,ied
        if (ISS%hmask(i,j) == 1) ISS%h_shelf(i,j) = h_after_flux2(i,j)
      enddo
    enddo

  endif

  if (CS%calc_flux_inout) call calculate_flux_inout(CS, ISS, G, uh_ice, vh_ice)

  if (CS%moving_shelf_front) then
    call shelf_advance_front(CS, ISS, G, ISS%hmask, uh_ice, vh_ice)
    if (CS%min_thickness_simple_calve > 0.0) then
      if (CS%use_DG_thickness) then
        call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                           CS%min_thickness_simple_calve, h_nodal=CS%h_nodal)
      else
        call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                           CS%min_thickness_simple_calve)
      endif
    endif
    if (CS%calve_to_mask) then
      if (CS%use_DG_thickness) then
        call calve_to_mask(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, CS%calve_mask, &
                           h_nodal=CS%h_nodal)
      else
        call calve_to_mask(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, CS%calve_mask)
      endif
    endif
  elseif (calve_ice_shelf_bergs) then
    !advect the front to create partially-filled cells
    call shelf_advance_front(CS, ISS, G, ISS%hmask, uh_ice, vh_ice, calving=calve_ice_shelf_bergs)
    !add mass of the partially-filled cells to calving field, which is used to initialize icebergs
    !Then, remove the partially-filled cells from the ice shelf
    ISS%calving(:,:) = 0.0
    ISS%calving_hflx(:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (ISS%hmask(i,j)==2) then
        ISS%calving(i,j) = (ISS%h_shelf(i,j) * CS%density_ice) * &
                           (ISS%area_shelf_h(i,j) * G%IareaT(i,j)) / time_step
        ISS%calving_hflx(i,j) = (CS%Cp_ice * CS%t_shelf(i,j)) * &
                                ((ISS%h_shelf(i,j) * CS%density_ice) * &
                                (ISS%area_shelf_h(i,j) * G%IareaT(i,j)))
        ISS%h_shelf(i,j) = 0.0 ; ISS%area_shelf_h(i,j) = 0.0 ; ISS%hmask(i,j) = 0.0
        if (CS%use_DG_thickness) then
          CS%h_nodal(i,j,:,:) = 0.0
        endif
      endif
    enddo ; enddo
  endif

  do j=jsc,jec ; do i=isc,iec
    ISS%mass_shelf(i,j) = ISS%h_shelf(i,j) * CS%density_ice
  enddo ; enddo

  call pass_var(ISS%mass_shelf, G%domain, complete=.false.)
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
  call pass_var(ISS%hmask, G%domain, complete=.true.)

  call update_velocity_masks(CS, G, ISS%hmask, CS%umask, CS%vmask, CS%u_face_mask, CS%v_face_mask)

end subroutine ice_shelf_advect

!>This subroutine computes u- and v-velocities of the ice shelf iterating on non-linear ice viscosity
!subroutine ice_shelf_solve_outer(CS, ISS, G, US, u_shlf, v_shlf, iters, time)
subroutine ice_shelf_solve_outer(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, iters, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  integer,                intent(out)   :: iters !< The number of iterations used in the solver.
  type(time_type),        intent(in)    :: Time !< The current model time

  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: taudx !< Driving x-stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: taudy !< Driving y-stress at q-points [R L3 Z T-2 ~> kg m s-2]
  !real, dimension(SZDIB_(G),SZDJB_(G)) :: u_bdry_cont ! Boundary u-stress contribution [R L3 Z T-2 ~> kg m s-2]
  !real, dimension(SZDIB_(G),SZDJB_(G)) :: v_bdry_cont ! Boundary v-stress contribution [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Au, Av ! The retarding lateral stress contributions [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: u_last, v_last ! Previous velocities [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: H_node ! Ice shelf thickness at corners [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Normvec  ! Velocities used for convergence [L2 T-2 ~> m2 s-2]
  logical :: converged ! Indicates nonlinear convergence
  logical :: calc_Au_for_convergence ! Used for convergence criteria than need a CG_action
  character(len=160) :: mesg  ! The text of an error message
  integer :: conv_flag, i, j, iter
  integer :: Isdq, Iedq, Jsdq, Jedq, isd, ied, jsd, jed
  integer :: Iscq, Iecq, Jscq, Jecq, isc, iec, jsc, jec
  real    :: err_max, err_tempu, err_tempv, err_init ! Errors in [R L3 Z T-2 ~> kg m s-2] or [L T-1 ~> m s-1]
  real    :: norm_tau, err_rr ! Errors in [R L3 Z T-2 ~> kg m s-2] for relative residual
  real    :: ew_resid       = 0.0  ! L2 norm of stress residual ||A(u)u - tau|| for Eisenstat-Walker [kg m s-2]
  real    :: ew_prev_resid  = 0.0  ! Previous ew_resid; 0.0 flags first Newton call [kg m s-2]
  real    :: ew_eta         = 0.0  ! Current EW inner tolerance [nondim]
  real    :: ew_eta_prev    = 0.0  ! Previous EW inner tolerance for Chacon 2008 sharp-decrease safeguard [nondim]
  real    :: ew_stol                ! Temporary safeguard tolerance [nondim]
  real    :: max_vel  ! The maximum velocity magnitude [L T-1 ~> m s-1]
  real    :: tempu, tempv   ! Temporary variables with velocity magnitudes [L T-1 ~> m s-1]
  real    :: Norm, PrevNorm ! Velocities used to assess convergence [L T-1 ~> m s-1]
  real    :: rhoi_rhow ! The density of ice divided by a typical water density [nondim]
  integer :: Is_sum, Js_sum, Ie_sum, Je_sum ! Loop bounds for global sums or arrays starting at 1.
  integer :: Iscq_sv, Jscq_sv ! Starting loop bound for sum_vec

  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  ! Determine the loop limits for sums, bearing in mind that the arrays will be starting at 1.
  ! Includes the edge of the tile is at the western/southern bdry (if symmetric)
  if (CS%nonlin_solve_err_mode >= 3 .or. CS%ssa_add_rel_resid) then
    if ((isc+G%idg_offset==G%isg) .and. (.not. CS%reentrant_x)) then
      Is_sum = Iscq + (1-Isdq) ; Iscq_sv = Iscq
    else
      Is_sum = isc  + (1-Isdq) ; Iscq_sv = isc
    endif
    if ((jsc+G%jdg_offset==G%jsg) .and. (.not. CS%reentrant_y)) then
      Js_sum = Jscq + (1-Jsdq) ; Jscq_sv = Jscq
    else
      Js_sum = jsc + (1-Jsdq) ; Jscq_sv = jsc
    endif
    Ie_sum = Iecq + (1-Isdq) ; Je_sum = Jecq + (1-Jsdq)
  endif

  taudx(:,:) = 0.0 ; taudy(:,:) = 0.0
  Au(:,:) = 0.0 ; Av(:,:) = 0.0

  ! need to make these conditional on GL interpolation
  H_node(:,:) = 0.0
  !CS%ground_frac(:,:) = 0.0

  if (.not. CS%GL_couple) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf) - CS%bed_elev(i,j) > 0) then
        CS%ground_frac(i,j) = 1.0
        CS%OD_av(i,j) =0.0
      endif
    enddo ; enddo
  endif

  ! Warning: This turns off Picard entirely and may not converge.
  if (CS%newton_after_tolerance<=0.0) CS%doing_newton=.true.

  ! Set CS%ground_frac in GL-regularize cells to the fraction of sub-grid integration
  ! points that are grounded (case 2: GL_regularize=True). Other cases leave ground_frac
  ! at the binary or running-mean value already set upstream. H_node is needed by the
  ! non-DG branch of compute_ground_frac and by CG_action_subgrid_basal further down.
  ! Computed before the driving-stress call so that the DG nsub switch and the non-DG
  ! Neumann test see the freshly-computed fractional ground_frac in the current outer
  ! iteration rather than lagged by one.
  if (CS%GL_regularize .and. .not. CS%use_DG_thickness) then
    call interpolate_H_to_B(G, ISS%h_shelf, ISS%hmask, H_node, CS%min_h_shelf)
  endif
  ! Refresh the continuous flotation-gate field before any grounded/floating
  ! decisions are made for this outer solve (h is frozen for its duration).
  if (CS%use_DG_thickness .and. CS%dg_gl_gate_continuous) call compute_h_flot(CS, ISS, G)
  call compute_ground_frac(CS, ISS, G, H_node)

  ! Calculate RHS
  if (CS%use_DG_thickness) then
    if (CS%dg_driving_stress_IBP) then
      call calc_shelf_driving_stress_DG(CS, ISS, G, US, taudx, taudy, CS%OD_av)
    else
      call calc_shelf_driving_stress_DG_strong(CS, ISS, G, US, taudx, taudy, CS%OD_av)
    endif
  else
    call calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, CS%OD_av)
  endif
  call pass_vector(taudx, taudy, G%domain, TO_ALL, BGRID_NE)

  ! Calculate basal drag constants and initial velocity
  call calc_shelf_basal_prefactors(CS, ISS, G, US)
  call calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
  if (CS%doing_newton) then
    call pass_var(CS%ice_visc, G%domain, complete=.false.)
    call pass_var(CS%newton_str_sh, G%domain, complete=.false.)
    call pass_var(CS%newton_visc_factor, G%domain, complete=.true.)
    call pass_vector(CS%newton_str_ux, CS%newton_str_vy, G%domain, TO_ALL, AGRID)
  else
    call pass_var(CS%ice_visc, G%domain, complete=.true.)
  endif

  ! Calculate err_init, the denominator for some convergence criteria
  if (CS%nonlin_solve_err_mode == 1 .or. CS%nonlin_solve_err_mode == 4) then
    Au(:,:) = 0.0 ; Av(:,:) = 0.0
    call CG_action(CS, Au, Av, u_shlf, v_shlf, CS%Phi, CS%Phisub, CS%umask, CS%vmask, ISS%hmask, H_node, &
      CS%ice_visc, CS%bed_elev, u_shlf, v_shlf, &
      G, US, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi_rhow, use_newton_in=.false., &
      h_shelf=ISS%h_shelf)
    call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE) ! TODO: is this needed?
  endif

  if (CS%nonlin_solve_err_mode == 1) then
    err_init = 0 ; err_tempu = 0 ; err_tempv = 0
    do J=G%JscB,G%JecB ; do I=G%IscB,G%IecB
      if (CS%umask(I,J) == 1) then
        err_tempu = ABS(Au(I,J) - taudx(I,J))
        if (err_tempu >= err_init) err_init = err_tempu
      endif
      if (CS%vmask(I,J) == 1) then
        err_tempv = ABS(Av(I,J) - taudy(I,J))
        if (err_tempv >= err_init) err_init = err_tempv
      endif
    enddo ; enddo
    call max_across_PEs(err_init)

  elseif (CS%nonlin_solve_err_mode == 3) then
    Normvec(:,:) = 0.0
    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
      if (CS%umask(I,J) == 1) Normvec(I,J) = (u_shlf(I,J)**2)
      if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + (v_shlf(I,J)**2)
    enddo ; enddo
    Norm = sqrt( reproducing_sum( Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, unscale=US%L_T_to_m_s**2 ) )

  elseif (CS%nonlin_solve_err_mode == 4) then
    Normvec(:,:) = 0.0
    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
      if (CS%umask(I,J) == 1) Normvec(I,J) = ((Au(I,J) - taudx(I,J))**2)
      if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + ((Av(I,J) - taudy(I,J))**2)
    enddo ; enddo
    err_init = sqrt(reproducing_sum(Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, &
      unscale=((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2))
  endif

  if (CS%nonlin_solve_err_mode == 5 .or. CS%ssa_add_rel_resid) then
    Normvec(:,:) = 0.0
    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
      if (CS%umask(I,J) == 1) Normvec(I,J) = (taudx(I,J)**2)
      if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + (taudy(I,J)**2)
    enddo ; enddo
    if (CS%nonlin_solve_err_mode == 5) then
      err_init = sqrt(reproducing_sum(Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, &
        unscale=((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2))
    else
      norm_tau = sqrt(reproducing_sum(Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, &
        unscale=((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2))
    endif
  endif

  u_last(:,:) = u_shlf(:,:) ; v_last(:,:) = v_shlf(:,:)
  if (CS%doing_newton) then
    CS%cg_tol_current = CS%cg_newton_tolerance
  else
    CS%cg_tol_current = CS%cg_tolerance
  endif
  ew_prev_resid  = 0.0
  converged = .false.
  calc_Au_for_convergence = (CS%nonlin_solve_err_mode == 1 .or. CS%nonlin_solve_err_mode == 4 .or. &
                             CS%nonlin_solve_err_mode == 5 .or. CS%ssa_add_rel_resid)

  !! begin loop

  do iter=1,50

    ! The linear solve
    call ice_shelf_solve_inner(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, H_node, &
                               ISS%hmask, conv_flag, iters, time, CS%Phi, CS%Phisub)

    if (CS%debug) then
      call qchksum(u_shlf, "u shelf", G%HI, haloshift=2, unscale=US%L_T_to_m_s)
      call qchksum(v_shlf, "v shelf", G%HI, haloshift=2, unscale=US%L_T_to_m_s)
    endif

    write(mesg,*) "ice_shelf_solve_outer: linear solve done in ",iters," iterations"
    call MOM_mesg(mesg, 5)

    ! Update viscosity
    call calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)

    if (CS%doing_newton) then
      call pass_var(CS%ice_visc, G%domain, complete=.false.)
      call pass_var(CS%newton_str_sh, G%domain, complete=.false.)
      call pass_var(CS%newton_visc_factor, G%domain, complete=.true.)
      call pass_vector(CS%newton_str_ux, CS%newton_str_vy, G%domain, TO_ALL, AGRID)
    else
      call pass_var(CS%ice_visc, G%domain, complete=.true.)
    endif

    ! Calculate convergence norms
    if (calc_Au_for_convergence) then
      Au(:,:) = 0 ; Av(:,:) = 0
      call CG_action(CS, Au, Av, u_shlf, v_shlf, CS%Phi, CS%Phisub, CS%umask, CS%vmask, ISS%hmask, &
        H_node, CS%ice_visc, CS%bed_elev, u_shlf, v_shlf, &
        G, US, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi_rhow, use_newton_in=.false., &
        h_shelf=ISS%h_shelf)

      if (CS%nonlin_solve_err_mode == 1) then
        err_max = 0

        do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
          if (CS%umask(I,J) == 1) then
            err_tempu = ABS(Au(I,J) - taudx(I,J))
            if (err_tempu >= err_max) err_max = err_tempu
          endif
          if (CS%vmask(I,J) == 1) then
            err_tempv = ABS(Av(I,J) - taudy(I,J))
            if (err_tempv >= err_max) err_max = err_tempv
          endif
        enddo ; enddo

        call max_across_PEs(err_max)
      endif

      if (CS%nonlin_solve_err_mode == 4 .or. CS%nonlin_solve_err_mode == 5 .or. CS%ssa_add_rel_resid) then
        Normvec(:,:) = 0.0
        do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
          if (CS%umask(I,J) == 1) Normvec(I,J) = ((Au(I,J) - taudx(I,J))**2)
          if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + ((Av(I,J) - taudy(I,J))**2)
        enddo ; enddo
        if (CS%nonlin_solve_err_mode == 4 .or. CS%nonlin_solve_err_mode == 5) then
          err_max = sqrt(reproducing_sum(Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, &
            unscale=((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2))
          if (CS%ssa_add_rel_resid) err_rr = err_max
        elseif (CS%ssa_add_rel_resid) then
          err_rr = sqrt(reproducing_sum(Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, &
            unscale=((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2))
        endif
      endif
    endif

    if (CS%nonlin_solve_err_mode == 2) then

      err_max=0. ;  max_vel = 0 ; tempu = 0 ; tempv = 0 ; err_tempu = 0
      do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
        if (CS%umask(I,J) == 1) then
          err_tempu = ABS(u_last(I,J)-u_shlf(I,J))
          if (err_tempu >= err_max) err_max = err_tempu
          tempu = u_shlf(I,J)
        else
          tempu = 0.0
        endif
        if (CS%vmask(I,J) == 1) then
          err_tempv = MAX(ABS(v_last(I,J)-v_shlf(I,J)), err_tempu)
          if (err_tempv >= err_max) err_max = err_tempv
          tempv = SQRT((v_shlf(I,J)**2) + (tempu**2))
        endif
        if (tempv >= max_vel) max_vel = tempv
      enddo ; enddo

      u_last(:,:) = u_shlf(:,:)
      v_last(:,:) = v_shlf(:,:)

      call max_across_PEs(max_vel)
      call max_across_PEs(err_max)
      err_init = max_vel

    elseif (CS%nonlin_solve_err_mode == 3) then
      PrevNorm = Norm ; Norm = 0.0 ; Normvec=0.0
      do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
        if (CS%umask(I,J) == 1) Normvec(I,J) = (u_shlf(I,J)**2)
        if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + (v_shlf(I,J)**2)
      enddo ; enddo
      Norm = sqrt( reproducing_sum( Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, unscale=US%L_T_to_m_s**2 ) )
      err_max = 2.*abs(Norm-PrevNorm) ; err_init = Norm+PrevNorm
    endif

    !Test convergence
    if (err_max <= CS%nonlinear_tolerance * err_init) then
      if (CS%ssa_add_rel_resid) then
        if (err_rr <= CS%rr_nonlinear_tolerance * norm_tau) converged = .true.
      else
        converged = .true.
      endif
    endif

    if (converged) then
      exit
    else
      write(mesg,*) "ice_shelf_solve_outer: nonlinear fractional residual = ", err_max/err_init
      call MOM_mesg(mesg, 5)

      if (CS%ssa_add_rel_resid) then
        write(mesg,*) "ice_shelf_solve_outer: nonlinear relative stress residual = ", err_rr/norm_tau
        call MOM_mesg(mesg, 5)
      endif

      ! Activate Newton
      if (err_max <= CS%newton_after_tolerance * err_init .and. .not. CS%doing_newton) then
        CS%doing_newton = .true.
        write(mesg,*) "ice_shelf_solve_outer: switching to Newton iterations at iter = ", iter
        call MOM_mesg(mesg, 7)
        call pass_var(CS%newton_str_sh, G%domain, complete=.false.)
        call pass_var(CS%newton_visc_factor, G%domain, complete=.true.)
        call pass_vector(CS%newton_str_ux, CS%newton_str_vy, G%domain, TO_ALL, AGRID)
        CS%cg_tol_current = CS%cg_newton_tolerance
      endif

      ! Inexact Newton: Adapt inner solver tolerance to prevent oversolving
      ! Based on Eisenstat-Walker Choice II (Eisenstat & Walker 1994): η_k = γ*(||F_k||/||F_{k-1}||)^α
      ! with γ=0.9, α=2 as default.  Uses the L2 norm of the nonlinear stress residual ||Au - tau||_2,
      ! consistent with the inner solver's convergence check (sv3dsums(3)).
      ! The first Newton step uses the standard cg_tolerance.
      if (CS%doing_newton .and. CS%newton_adapt_cg_tol) then
        !calculate residual needed for EW; some convergence criteria already did this
        if (CS%nonlin_solve_err_mode >= 4) then
          ew_resid=err_max
        elseif (CS%ssa_add_rel_resid) then
          ew_resid=err_rr
        else
          if (.not. calc_Au_for_convergence) then
            Au(:,:) = 0 ; Av(:,:) = 0
            call CG_action(CS, Au, Av, u_shlf, v_shlf, CS%Phi, CS%Phisub, CS%umask, CS%vmask, ISS%hmask, &
              H_node, CS%ice_visc, CS%bed_elev, u_shlf, v_shlf, &
              G, US, G%isc-1, G%iec+1, G%jsc-1, G%jec+1, rhoi_rhow, use_newton_in=.false., &
              h_shelf=ISS%h_shelf)
          endif
          Normvec(:,:) = 0.0
          do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
            if (CS%umask(I,J) == 1) Normvec(I,J) = ((Au(I,J) - taudx(I,J))**2)
            if (CS%vmask(I,J) == 1) Normvec(I,J) = Normvec(I,J) + ((Av(I,J) - taudy(I,J))**2)
          enddo ; enddo
          ew_resid = sqrt(reproducing_sum(Normvec, Is_sum, Ie_sum, Js_sum, Je_sum, &
            unscale=((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2))
        endif

        if (ew_prev_resid == 0.0) then
          ! First Newton iteration: seed residuals; use initial newton cg_tolerance this step
          ew_prev_resid  = ew_resid
          CS%cg_tol_current = CS%cg_newton_tolerance
          ew_eta_prev = CS%cg_tol_current
        else
          ! Safeguarding and oversolving adjustments:
          ! Eisenstat-Walker Choice II safeguard base formula
          ew_eta = CS%ew_gamma * (ew_resid / ew_prev_resid)**CS%ew_alpha
          ew_stol = CS%ew_gamma * ew_eta_prev**CS%ew_alpha
          !Safeguards to sharp decrease/oversolving:
          if (CS%ew_safety==1) then
            ! Eisenstat-Walker Choice II safeguard:
            write(mesg,*) "ice_shelf_solve_outer: ew_stol = ", ew_stol
            call MOM_mesg(mesg, 8)
            if (ew_stol > CS%ew_1_thres) ew_eta = max(ew_eta, ew_stol)
          elseif (CS%ew_safety==2) then
            ! PETSc choice 3 safeguard (e,g, Chacon 2008, J. Phys: Conf. Ser. 125 012041):
            ! Avoid steep decreases in ew_eta
            ew_eta  = min(CS%cg_newton_tolerance, max(ew_eta, ew_stol))
            ! Avoid oversolving in last Newton iters:
            ! The original is technically only applicable for nonlin_solve_err_mode=4:
            ! ew_stol = CS%ew_gamma * ew_resid_first * CS%nonlinear_tolerance / ew_resid
            ! Here, adapt for all nonlin_solve_err_modes:
            ew_stol = CS%ew_gamma * err_init * CS%nonlinear_tolerance / err_max
            if (CS%ssa_add_rel_resid) then
              ew_stol = min(ew_stol, CS%ew_gamma * norm_tau * CS%rr_nonlinear_tolerance / err_rr)
            endif
            ew_eta  = min(CS%cg_newton_tolerance, max(ew_eta, ew_stol))
            write(mesg,*) "ice_shelf_solve_outer: ew_stol = ", ew_stol
            call MOM_mesg(mesg, 8)
          endif
          ew_eta = min(ew_eta,CS%ew_eta_max)
          CS%cg_tol_current = ew_eta
          ew_eta_prev   = ew_eta
          ew_prev_resid = ew_resid
          write(mesg,*) "ice_shelf_solve_outer: New inner tolerance = ", CS%cg_tol_current
          call MOM_mesg(mesg, 8)
        endif
      endif
    endif
  enddo
  CS%doing_newton = .false.
  CS%cg_tol_current = CS%cg_tolerance

  write(mesg,*) "ice_shelf_solve_outer: nonlinear fractional residual = ", err_max/err_init
  call MOM_mesg(mesg)
  if (CS%ssa_add_rel_resid) then
    write(mesg,*) "ice_shelf_solve_outer: nonlinear relative residual = ", err_rr/norm_tau
    call MOM_mesg(mesg, 5)
  endif
  write(mesg,*) "ice_shelf_solve_outer: exiting nonlinear solve after ",iter," iterations"
  call MOM_mesg(mesg)

end subroutine ice_shelf_solve_outer

!> Unified inner linear solver for ice shelf velocity.
!! Performs shared setup (RHS, preconditioner, initial matrix-vector product),
!! dispatches to the selected Krylov method, and applies boundary conditions.
subroutine ice_shelf_solve_inner(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, H_node, &
                                  hmask, conv_flag, iters, time, Phi, Phisub)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: taudx !< The x-direction driving stress [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: taudy  !< The y-direction driving stress [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal (corner) points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                                 !! partly or fully covered by an ice-shelf
  integer,                intent(out)   :: conv_flag !< A flag indicating whether (1) or not (0) the
                                                     !! iterations have converged to the specified tolerance
  integer,                intent(out)   :: iters !< The number of iterations used in the solver.
  type(time_type),        intent(in)    :: Time !< The current model time
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< The gradients of bilinear basis elements at Gaussian
                                               !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                                  !! locations for finite element calculations [nondim]

  real, dimension(SZDIB_(G),SZDJB_(G)) :: &
        RHSu, RHSv, &      ! Right hand side of the stress balance [R L3 Z T-2 ~> m kg s-2]
        Au, Av, &          ! Matrix-vector product A*x [R L3 Z T-2 ~> kg m s-2]
        DIAGu, DIAGv, &    ! Diagonals [R L2 Z T-1 ~> kg s-1]
        IDIAGu, IDIAGv     ! Reciprocal diagonals [R-1 L-2 Z-1 T ~> kg-1 s]
  real    :: rhoi_rhow     ! The density of ice divided by a typical water density [nondim]
  real    :: resid_scale   ! A scaling factor for redimensionalizing the global residuals
                           ! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]
  integer :: Is_sum, Js_sum, Ie_sum, Je_sum ! Loop bounds for global sums or arrays starting at 1.
  integer :: Iscq_sv, Jscq_sv ! Starting loop bound for sum_vec arrays
  integer :: I, J
  integer :: Isdq, Iedq, Jsdq, Jedq, Iscq, Iecq, Jscq, Jecq
  integer :: isc, iec, jsc, jec

  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  ! Initialize shared arrays
  Au(:,:) = 0 ; Av(:,:) = 0 ; DIAGu(:,:) = 0 ; DIAGv(:,:) = 0

  ! Determine the loop limits for sums, bearing in mind that the arrays will be starting at 1.
  ! Includes the edge of the tile is at the western/southern bdry (if symmetric)
  if ((isc+G%idg_offset==G%isg) .and. (.not. CS%reentrant_x)) then
    Is_sum = Iscq + (1-Isdq) ; Iscq_sv = Iscq
  else
    Is_sum = isc  + (1-Isdq) ; Iscq_sv = isc
  endif
  if ((jsc+G%jdg_offset==G%jsg) .and. (.not. CS%reentrant_y)) then
    Js_sum = Jscq + (1-Jsdq) ; Jscq_sv = Jscq
  else
    Js_sum = jsc + (1-Jsdq) ; Jscq_sv = jsc
  endif
  Ie_sum = Iecq + (1-Isdq) ; Je_sum = Jecq + (1-Jsdq)

  RHSu(:,:) = taudx(:,:) ; RHSv(:,:) = taudy(:,:)
  call pass_vector(RHSu, RHSv, G%domain, TO_ALL, BGRID_NE, complete=.false.)

  call matrix_diagonal(CS, G, US, H_node, CS%ice_visc, u_shlf, v_shlf, &
                       hmask, rhoi_rhow, Phi, Phisub, DIAGu, DIAGv, h_shelf=ISS%h_shelf)
  call pass_vector(DIAGu, DIAGv, G%domain, TO_ALL, BGRID_NE, complete=.false.)

  call CG_action(CS, Au, Av, u_shlf, v_shlf, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                 H_node, CS%ice_visc, CS%bed_elev, u_shlf, v_shlf, &
                 G, US, isc-1, iec+1, jsc-1, jec+1, rhoi_rhow, use_newton_in=.false., &
                 h_shelf=ISS%h_shelf)
  call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE, complete=.true.)

  ! Precompute reciprocal diagonal
  IDIAGu(:,:) = 0.0 ; IDIAGv(:,:) = 0.0
  do J=Jsdq,Jedq ; do I=Isdq,Iedq
    if (CS%umask(I,J)==1 .AND. DIAGu(I,J)/=0) IDIAGu(I,J) = 1.0 / DIAGu(I,J)
    if (CS%vmask(I,J)==1 .AND. DIAGv(I,J)/=0) IDIAGv(I,J) = 1.0 / DIAGv(I,J)
  enddo ; enddo

  resid_scale = US%s_to_T*(US%RZL2_to_kg*US%L_T_to_m_s**2)

  ! Dispatch to selected solver
  select case (CS%inner_solver)
    case (INNER_CG)
      call ice_shelf_solve_inner_CG(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                    IDIAGu, IDIAGv, H_node, hmask, &
                                    rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                    Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv, &
                                    h_shelf=ISS%h_shelf)
    case (INNER_MINRES)
      call ice_shelf_solve_inner_MINRES(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                        IDIAGu, IDIAGv, H_node, hmask, &
                                        rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                        Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv, &
                                        h_shelf=ISS%h_shelf)
    case (INNER_CR)
      call ice_shelf_solve_inner_CR(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                    IDIAGu, IDIAGv, H_node, hmask, &
                                    rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                    Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv, &
                                    h_shelf=ISS%h_shelf)
  end select

  ! Apply boundary conditions
  do J=Jsdq,Jedq ; do I=Isdq,Iedq
      if (CS%umask(I,J) == 3) then
        u_shlf(I,J) = CS%u_bdry_val(I,J)
      elseif (CS%umask(I,J) == 0) then
        u_shlf(I,J) = 0
      endif

      if (CS%vmask(I,J) == 3) then
        v_shlf(I,J) = CS%v_bdry_val(I,J)
      elseif (CS%vmask(I,J) == 0) then
        v_shlf(I,J) = 0
      endif
  enddo ; enddo

  call pass_vector(u_shlf, v_shlf, G%domain, TO_ALL, BGRID_NE)

  if (conv_flag == 0) then
    iters = CS%cg_max_iterations
  endif

end subroutine ice_shelf_solve_inner

!> CG (Conjugate Gradient) inner Krylov solve for ice shelf velocity.
subroutine ice_shelf_solve_inner_CG(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                     IDIAGu, IDIAGv, H_node, hmask, &
                                     rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                     Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv, &
                                     h_shelf)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSu !< Right hand side, x [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSv !< Right hand side, y [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Au !< Matrix-vector product workspace, x [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Av !< Matrix-vector product workspace, y [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGu !< Reciprocal Jacobi diagonal, x [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGv !< Reciprocal Jacobi diagonal, y [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal points [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf coverage mask
  real,                   intent(in)    :: rhoi_rhow !< Ice-to-ocean density ratio [nondim]
  real,                   intent(in)    :: resid_scale !< Scaling for inner products
                                                       !! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< Basis element gradients at quadrature points [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Subgridscale quadrature weights [nondim]
  integer,                intent(out)   :: conv_flag !< Convergence flag: 1=converged, 0=not
  integer,                intent(out)   :: iters !< The number of iterations used
  integer,                intent(in)    :: Is_sum !< Starting i-index for global sums
  integer,                intent(in)    :: Js_sum !< Starting j-index for global sums
  integer,                intent(in)    :: Ie_sum !< Ending i-index for global sums
  integer,                intent(in)    :: Je_sum !< Ending j-index for global sums
  integer,                intent(in)    :: Iscq_sv !< Starting i-index for sum_vec arrays
  integer,                intent(in)    :: Jscq_sv !< Starting j-index for sum_vec arrays
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          optional, intent(in) :: h_shelf !< Ice shelf thickness on tracer grid [Z ~> m]

  real, dimension(SZDIB_(G),SZDJB_(G))  :: u_curr  !< Frozen current iterate u^k, used to evaluate basal friction
                                                   !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G))  :: v_curr  !< Frozen current iterate v^k, used to evaluate basal friction
                                                   !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) ::  &
                        Ru, Rv, &     ! Residuals [R L3 Z T-2 ~> m kg s-2]
                        Zu, Zv, &     ! Preconditioned residuals [L T-1 ~> m s-1]
                        Du, Dv        ! Search directions [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: sum_vec ! Pointwise D·A products for the alpha_k global sum
                                                  ! [kg m2 s-3]
  real, dimension(SZDIB_(G),SZDJB_(G),2) :: sum_vec_3d ! Array used for various residuals
                                                       ! sum_vec_3d(:,:,1) [kg m2 s-3]
                                                       ! sum_vec_3d(:,:,2) [kg2 m2 s-4]
  real    :: beta_k      ! Ratio of residuals used to update search direction [nondim]
  real    :: resid0tol2  ! Convergence tolerance times the initial residual [m2 kg2 s-4]
  real    :: sv3dsum     ! An unused variable returned when taking global sum of residuals [various]
  real    :: sv3dsums(2) ! The index-wise global sums of sum_vec_3d
                         ! sv3dsums(1) [kg m2 s-3]
                         ! sv3dsums(2) [kg2 m2 s-4]
  real    :: alpha_k     ! A scaling factor for iterative corrections [nondim]
  real    :: rho_old     ! The preconditioned residual inner product Z·R from the previous CG
                         ! iteration, scaled by resid_scale [kg m2 s-3]
  real    :: resid2_scale ! A scaling factor for redimensionalizing the global squared residuals
                          ! [T4 kg2 m2 R-2 Z-2 L-6 s-4 ~> 1]
  integer :: cg_halo     ! Number of halo vertices to include during a CG iteration
  integer :: max_cg_halo ! Maximum possible number of halo vertices to include in the CG iterations
  integer :: iter, i, j, isc, iec, jsc, jec, is, js, ie, je, is2, ie2, js2, je2
  integer :: Isdq, Iedq, Jsdq, Jedq, Iscq, Iecq, Jscq, Jecq, nx_halo, ny_halo

  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  ny_halo = G%domain%njhalo ; nx_halo = G%domain%nihalo
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  resid2_scale = ((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2

  Ru(:,:) = 0 ; Rv(:,:) = 0 ; Zu(:,:) = 0 ; Zv(:,:) = 0 ; Du(:,:) = 0 ; Dv(:,:) = 0

  Ru(:,:) = (RHSu(:,:) - Au(:,:)) ; Rv(:,:) = (RHSv(:,:) - Av(:,:))

  ! current velocities used in CG_action for basal drag
  u_curr(:,:) = u_shlf(:,:) ; v_curr(:,:) = v_shlf(:,:)

  do J=Jsdq,Jedq ; do I=Isdq,Iedq
    if (CS%umask(I,J) == 1) Zu(I,J) = Ru(I,J) * IDIAGu(I,J)
    if (CS%vmask(I,J) == 1) Zv(I,J) = Rv(I,J) * IDIAGv(I,J)
    Du(I,J) = Zu(I,J)
    Dv(I,J) = Zv(I,J)
  enddo ; enddo

  ! Compute rho_old = Z·R and resid0tol2 before the CG loop
  sum_vec_3d(:,:,:) = 0.0
  do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
    if (CS%umask(I,J) == 1) then
      sum_vec_3d(I,J,1) = resid_scale  * (Zu(I,J) * Ru(I,J))
      sum_vec_3d(I,J,2) = resid2_scale * Ru(I,J)**2
    endif
    if (CS%vmask(I,J) == 1) then
      sum_vec_3d(I,J,1) = sum_vec_3d(I,J,1) + resid_scale  * (Zv(I,J) * Rv(I,J))
      sum_vec_3d(I,J,2) = sum_vec_3d(I,J,2) + resid2_scale * Rv(I,J)**2
    endif
  enddo ; enddo

  sv3dsum = reproducing_sum( sum_vec_3d(:,:,1:2), Is_sum, Ie_sum, Js_sum, Je_sum, sums=sv3dsums(1:2) )

  rho_old = sv3dsums(1)
  !resid0 = sqrt(sv3dsums(2))
  resid0tol2 = CS%cg_tol_current**2 * sv3dsums(2)

  if (G%symmetric) then
    max_cg_halo=min(nx_halo,ny_halo)
  else
    max_cg_halo=min(nx_halo,ny_halo)-1
  endif
  cg_halo = max_cg_halo
  conv_flag = 0

  if (CS%cg_halo_shrink) then
    is = isc - cg_halo ; ie = Iecq + cg_halo
    js = jsc - cg_halo ; je = Jecq + cg_halo
    is2 = is ; ie2 = ie-1
    js2 = js ; je2 = je-1
  else
    is = isc - 1 ; ie = iec + 1
    js = jsc - 1 ; je = jec + 1
    is2 = Iscq ; ie2 = Iecq
    js2 = Jscq ; je2 = Jecq
  endif

  !!!!!!!!!!!!!!!!!!
  !!              !!
  !! MAIN CG LOOP !!
  !!              !!
  !!!!!!!!!!!!!!!!!!

  do iter = 1,CS%cg_max_iterations

    Au(:,:) = 0 ; Av(:,:) = 0

    call CG_action(CS, Au, Av, Du, Dv, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                   H_node, CS%ice_visc, CS%bed_elev, u_curr, v_curr, &
                   G, US, is, ie, js, je, rhoi_rhow, h_shelf=h_shelf)

    sum_vec(:,:) = 0.0

    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
      if (CS%umask(I,J) == 1) sum_vec(I,J) = resid_scale * (Du(I,J) * Au(I,J))
      if (CS%vmask(I,J) == 1) sum_vec(I,J) = sum_vec(I,J) + resid_scale * (Dv(I,J) * Av(I,J))
    enddo ; enddo

    sv3dsum = reproducing_sum( sum_vec(:,:), Is_sum, Ie_sum, Js_sum, Je_sum )

    if (sv3dsum == 0.0) then
      iters = iter
      conv_flag = 1
      exit
    endif

    alpha_k = rho_old / sv3dsum

    do J=js2,je2 ; do I=is2,ie2
      if (CS%umask(I,J) == 1) then
        u_shlf(I,J) = u_shlf(I,J) + alpha_k * Du(I,J)
        Ru(I,J) = Ru(I,J) - alpha_k * Au(I,J)
        Zu(I,J) = Ru(I,J) * IDIAGu(I,J)
      endif
      if (CS%vmask(I,J) == 1) then
        v_shlf(I,J) = v_shlf(I,J) + alpha_k * Dv(I,J)
        Rv(I,J) = Rv(I,J) - alpha_k * Av(I,J)
        Zv(I,J) = Rv(I,J) * IDIAGv(I,J)
      endif
    enddo ; enddo

    ! beta_k = (Z \dot R) / (Z_prev \dot R_prev)
    sum_vec_3d(:,:,:) = 0.0 ; sv3dsums(:)=0.0

    do J=jscq_sv,jecq ; do i=iscq_sv,iecq
      if (CS%umask(I,J) == 1) then
        sum_vec_3d(I,J,1) = resid_scale  * (Zu(I,J) * Ru(I,J))
        sum_vec_3d(I,J,2) = resid2_scale * Ru(I,J)**2
      endif
      if (CS%vmask(I,J) == 1) then
        sum_vec_3d(I,J,1) = sum_vec_3d(I,J,1) + resid_scale  * (Zv(I,J) * Rv(I,J))
        sum_vec_3d(I,J,2) = sum_vec_3d(I,J,2) + resid2_scale * Rv(I,J)**2
      endif
    enddo ; enddo

    sv3dsum = reproducing_sum( sum_vec_3d(:,:,1:2), Is_sum, Ie_sum, Js_sum, Je_sum, sums=sv3dsums(1:2) )

    beta_k = sv3dsums(1) / rho_old

    if (sv3dsums(2) <= resid0tol2) then
      iters = iter
      conv_flag = 1
      exit
    endif

    do J=js2,je2 ; do I=is2,ie2
      if (CS%umask(I,J) == 1) Du(I,J) = Zu(I,J) + beta_k * Du(I,J)
      if (CS%vmask(I,J) == 1) Dv(I,J) = Zv(I,J) + beta_k * Dv(I,J)
    enddo ; enddo

    rho_old = sv3dsums(1)

    if (CS%cg_halo_shrink) then
      cg_halo = cg_halo - 1
      if (cg_halo == 0) then
        call pass_vector(Du, Dv, G%domain, TO_ALL, BGRID_NE, complete=.false.)
        call pass_vector(Zu, Zv, G%domain, TO_ALL, BGRID_NE, complete=.false.)
        call pass_vector(Ru, Rv, G%domain, TO_ALL, BGRID_NE, complete=.false.)
        call pass_vector(u_shlf, v_shlf, G%domain, TO_ALL, BGRID_NE, complete=.true.)
        cg_halo = max_cg_halo
      endif
      is = isc - cg_halo ; ie = Iecq + cg_halo
      js = jsc - cg_halo ; je = Jecq + cg_halo
      is2 = is ; ie2 = ie-1
      js2 = js ; je2 = je-1
    else
      call pass_vector(Du, Dv, G%domain, TO_ALL, BGRID_NE)
    endif

  enddo ! end of CG loop

end subroutine ice_shelf_solve_inner_CG

!> MINRES inner Krylov solve for ice shelf velocity.
subroutine ice_shelf_solve_inner_MINRES(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                         IDIAGu, IDIAGv, H_node, hmask, &
                                         rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                         Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv, &
                                         h_shelf)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSu !< Right hand side, x [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSv !< Right hand side, y [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Au !< Matrix-vector product workspace, x [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Av !< Matrix-vector product workspace, y [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGu !< Reciprocal Jacobi diagonal, x [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGv !< Reciprocal Jacobi diagonal, y [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal points [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf coverage mask
  real,                   intent(in)    :: rhoi_rhow !< Ice-to-ocean density ratio [nondim]
  real,                   intent(in)    :: resid_scale !< Scaling for inner products
                                                       !! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< Basis element gradients at quadrature points [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Subgridscale quadrature weights [nondim]
  integer,                intent(out)   :: conv_flag !< Convergence flag: 1=converged, 0=not
  integer,                intent(out)   :: iters !< The number of iterations used
  integer,                intent(in)    :: Is_sum !< Starting i-index for global sums
  integer,                intent(in)    :: Js_sum !< Starting j-index for global sums
  integer,                intent(in)    :: Ie_sum !< Ending i-index for global sums
  integer,                intent(in)    :: Je_sum !< Ending j-index for global sums
  integer,                intent(in)    :: Iscq_sv !< Starting i-index for sum_vec arrays
  integer,                intent(in)    :: Jscq_sv !< Starting j-index for sum_vec arrays
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          optional, intent(in) :: h_shelf !< Ice shelf thickness on tracer grid [Z ~> m]

  real, dimension(SZDIB_(G),SZDJB_(G))  :: u_curr  !< Frozen current iterate u^k, used to evaluate basal friction
                                                   !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G))  :: v_curr  !< Frozen current iterate v^k, used to evaluate basal friction
                                                   !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) ::  &
        V_old_u, V_old_v, V_curr_u, V_curr_v, V_new_u, V_new_v, & ! Lanczos basis vectors [R L3 Z T-2 ~> m kg s-2]
        Z_curr_u, Z_curr_v, Z_new_u, Z_new_v, &    ! Preconditioned Lanczos vectors [L T-1 ~> m s-1]
        W_old_u, W_old_v, W_curr_u, W_curr_v, W_new_u, W_new_v, & ! MINRES search directions [L T-1 ~> m s-1]
        Qu, Qv            ! A * Z_curr [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: sum_vec_3d ! Pointwise products for global sums
                                                     ! [kg m2 s-3] before normalization;
                                                     ! [nondim] inside loop (after Lanczos normalization)
  real    :: alpha       ! Lanczos diagonal element (Rayleigh quotient) [nondim]
  real    :: beta1       ! Current Lanczos off-diagonal coefficient;
                         ! initial value [kg^1/2 m s^-3/2], then [nondim] after iter 1
  real    :: beta2       ! Next Lanczos off-diagonal coefficient [nondim]
  real    :: eta         ! MINRES residual norm estimate [kg^1/2 m s^-3/2]
  real    :: eta_curr    ! Effective step magnitude for current iteration [kg^1/2 m s^-3/2]
  real    :: c0, s0, c1, s1, c2, s2  ! Givens rotation cosines and sines [nondim]
  real    :: d0, d1, d2  ! Tridiagonal QR factorization coefficients [nondim]
  real    :: resid0tol   ! Convergence tolerance (CS%cg_tol_newton * beta1) [kg^1/2 m s^-3/2]
  real    :: current_norm ! Current MINRES residual norm estimate [kg^1/2 m s^-3/2]
  real    :: sv3dsum     ! Global reproducing sum of sum_vec_3d;
                         ! [kg m2 s-3] before normalization, [nondim] inside loop
  real    :: Ibeta1      ! Reciprocal of initial beta1 [kg^-1/2 m-1 s^3/2]
  real    :: Ibeta2      ! Reciprocal of beta2 [nondim]
  real    :: Id1         ! Reciprocal of d1 [nondim]
  integer :: iter, i, j, isc, iec, jsc, jec
  integer :: Isdq, Iedq, Jsdq, Jedq, Iscq, Iecq, Jscq, Jecq

  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  ! Initialize MINRES-specific arrays
  V_old_u(:,:) = 0 ; V_old_v(:,:) = 0 ; V_curr_u(:,:) = 0 ; V_curr_v(:,:) = 0
  Z_curr_u(:,:) = 0 ; Z_curr_v(:,:) = 0
  W_old_u(:,:) = 0 ; W_old_v(:,:) = 0 ; W_curr_u(:,:) = 0 ; W_curr_v(:,:) = 0
  Qu(:,:) = 0 ; Qv(:,:) = 0

  ! Initial Residual
  V_curr_u(:,:) = (RHSu(:,:) - Au(:,:)) ; V_curr_v(:,:) = (RHSv(:,:) - Av(:,:))

  ! current velocities used in CG_action for basal drag
  u_curr(:,:) = u_shlf(:,:) ; v_curr(:,:) = v_shlf(:,:)

  do J=Jscq,Jecq ; do I=Iscq,Iecq
     if (CS%umask(I,J) == 1) Z_curr_u(I,J) = V_curr_u(I,J) * IDIAGu(I,J)
     if (CS%vmask(I,J) == 1) Z_curr_v(I,J) = V_curr_v(I,J) * IDIAGv(I,J)
  enddo ; enddo

  sum_vec_3d(:,:) = 0.0
  do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
    if (CS%umask(I,J) == 1) sum_vec_3d(I,J) = resid_scale * (V_curr_u(I,J) * Z_curr_u(I,J))
    if (CS%vmask(I,J) == 1) sum_vec_3d(I,J) = sum_vec_3d(I,J) + resid_scale * (V_curr_v(I,J) * Z_curr_v(I,J))
  enddo ; enddo
  sv3dsum = reproducing_sum( sum_vec_3d(:,:), Is_sum, Ie_sum, Js_sum, Je_sum )

  beta1 = sqrt(abs(sv3dsum))

  if (beta1 == 0.0) then
     conv_flag = 1
     iters = 0
     return
  endif

  Ibeta1 = 1.0/beta1

  ! Normalize initial Lanczos vectors
  do J=Jscq,Jecq ; do I=Iscq,Iecq
     if (CS%umask(I,J) == 1) then
         V_curr_u(I,J) = V_curr_u(I,J) * Ibeta1
         Z_curr_u(I,J) = Z_curr_u(I,J) * Ibeta1
     endif
     if (CS%vmask(I,J) == 1) then
         V_curr_v(I,J) = V_curr_v(I,J) * Ibeta1
         Z_curr_v(I,J) = Z_curr_v(I,J) * Ibeta1
     endif
  enddo ; enddo

  ! Sync Z_curr prior to entering the loop
  call pass_vector(Z_curr_u, Z_curr_v, G%domain, TO_ALL, BGRID_NE)

  eta = beta1
  resid0tol = CS%cg_tol_current * beta1
  conv_flag = 0

  c0 = 1.0 ; s0 = 0.0 ; c1 = 1.0 ; s1 = 0.0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                !!
  !! MAIN MINRES LANCZOS LOOP       !!
  !!                                !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do iter = 1, CS%cg_max_iterations

    ! --- STEP 1: Matrix Vector Product ---
    Qu(:,:) = 0 ; Qv(:,:) = 0
    call CG_action(CS, Qu, Qv, Z_curr_u, Z_curr_v, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                   H_node, CS%ice_visc, CS%bed_elev, u_curr, v_curr, &
                   G, US, isc-1, iec+1, jsc-1, jec+1, rhoi_rhow, h_shelf=h_shelf)
    ! --- STEP 2: alpha = q dot z_curr ---
    sum_vec_3d(:,:) = 0.0
    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
       if (CS%umask(I,J) == 1) sum_vec_3d(I,J) = resid_scale * (Qu(I,J) * Z_curr_u(I,J))
       if (CS%vmask(I,J) == 1) sum_vec_3d(I,J) = sum_vec_3d(I,J) + resid_scale * (Qv(I,J) * Z_curr_v(I,J))
    enddo ; enddo
    sv3dsum = reproducing_sum( sum_vec_3d(:,:), Is_sum, Ie_sum, Js_sum, Je_sum )
    alpha = sv3dsum

    ! --- FUSED STEPS 3 & 4: Update V_new and Precondition to Z_new ---
    do J=Jscq,Jecq ; do I=Iscq,Iecq
       if (CS%umask(I,J) == 1) then
           V_new_u(I,J) = Qu(I,J) - alpha * V_curr_u(I,J) - beta1 * V_old_u(I,J)
           Z_new_u(I,J) = V_new_u(I,J) * IDIAGu(I,J)
       endif
       if (CS%vmask(I,J) == 1) then
           V_new_v(I,J) = Qv(I,J) - alpha * V_curr_v(I,J) - beta1 * V_old_v(I,J)
           Z_new_v(I,J) = V_new_v(I,J) * IDIAGv(I,J)
       endif
    enddo ; enddo

    ! --- STEP 5: beta2 = sqrt(v_new dot z_new) ---
    sum_vec_3d(:,:) = 0.0
    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
       if (CS%umask(I,J) == 1) sum_vec_3d(I,J) = resid_scale * (V_new_u(I,J) * Z_new_u(I,J))
       if (CS%vmask(I,J) == 1) sum_vec_3d(I,J) = sum_vec_3d(I,J) + resid_scale * (V_new_v(I,J) * Z_new_v(I,J))
    enddo ; enddo
    sv3dsum = reproducing_sum( sum_vec_3d(:,:), Is_sum, Ie_sum, Js_sum, Je_sum )
    beta2 = sqrt(abs(sv3dsum))

    ! --- STEP 6: Apply Givens Rotations ---
    d0 = c1 * alpha - c0 * s1 * beta1
    d1 = sqrt(d0**2 + beta2**2)

    if (d1 == 0.0) then
      iters = iter
      conv_flag = 1
      exit
    endif

    Id1 = 1.0 / d1
    if (beta2 > 0) Ibeta2 = 1.0 / beta2

    d2 = s1 * alpha + c0 * c1 * beta1
    c2 = d0 * Id1
    s2 = beta2 * Id1

    eta_curr = c2 * eta
    eta = -s2 * eta
    current_norm = abs(eta)

    ! --- FUSED STEPS 7 & 9: Update u/v, Check Convergence, and Shift Vectors ---
    do J=Jscq,Jecq ; do I=Iscq,Iecq
       if (CS%umask(I,J) == 1) then
           W_new_u(I,J) = (Z_curr_u(I,J) - (d2 * W_curr_u(I,J) + beta1 * s0 * W_old_u(I,J))) * Id1
           u_shlf(I,J) = u_shlf(I,J) + eta_curr * W_new_u(I,J)
           if (beta2 > 0.0) then
               V_old_u(I,J) = V_curr_u(I,J)
               V_curr_u(I,J) = V_new_u(I,J) * Ibeta2
               Z_curr_u(I,J) = Z_new_u(I,J) * Ibeta2
               W_old_u(I,J) = W_curr_u(I,J)
               W_curr_u(I,J) = W_new_u(I,J)
           endif
       endif
       if (CS%vmask(I,J) == 1) then
           W_new_v(I,J) = (Z_curr_v(I,J) - (d2 * W_curr_v(I,J) + beta1 * s0 * W_old_v(I,J))) * Id1
           v_shlf(I,J) = v_shlf(I,J) + eta_curr * W_new_v(I,J)
           if (beta2 > 0.0) then
               V_old_v(I,J) = V_curr_v(I,J)
               V_curr_v(I,J) = V_new_v(I,J) * Ibeta2
               Z_curr_v(I,J) = Z_new_v(I,J) * Ibeta2
               W_old_v(I,J) = W_curr_v(I,J)
               W_curr_v(I,J) = W_new_v(I,J)
           endif
       endif
    enddo ; enddo

    ! --- STEP 8: Check Convergence ---
    if (current_norm <= resid0tol .or. beta2 == 0.0) then
      iters = iter
      conv_flag = 1
      exit
    endif

    ! Sync Z_curr for the next iteration's CG_action
    call pass_vector(Z_curr_u, Z_curr_v, G%domain, TO_ALL, BGRID_NE)

    beta1 = beta2
    c0 = c1 ; c1 = c2
    s0 = s1 ; s1 = s2

  enddo ! end of MINRES loop

end subroutine ice_shelf_solve_inner_MINRES

!> CR (Conjugate Residual) inner Krylov solve for ice shelf velocity.
subroutine ice_shelf_solve_inner_CR(CS, G, US, u_shlf, v_shlf, RHSu, RHSv, Au, Av, &
                                     IDIAGu, IDIAGv, H_node, hmask, &
                                     rhoi_rhow, resid_scale, Phi, Phisub, conv_flag, iters, &
                                     Is_sum, Js_sum, Ie_sum, Je_sum, Iscq_sv, Jscq_sv, &
                                     h_shelf)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_shlf  !< The zonal ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_shlf  !< The meridional ice shelf velocity [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSu !< Right hand side, x [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: RHSv !< Right hand side, y [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Au !< Matrix-vector product workspace, x [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: Av !< Matrix-vector product workspace, y [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGu !< Reciprocal Jacobi diagonal, x [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: IDIAGv !< Reciprocal Jacobi diagonal, y [R-1 L-2 Z-1 T ~> kg-1 s]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal points [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf coverage mask
  real,                   intent(in)    :: rhoi_rhow !< Ice-to-ocean density ratio [nondim]
  real,                   intent(in)    :: resid_scale !< Scaling for inner products
                                                       !! [T3 kg m2 R-1 Z-1 L-4 s-3 ~> 1]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< Basis element gradients at quadrature points [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), &
                          intent(in)    :: Phisub !< Subgridscale quadrature weights [nondim]
  integer,                intent(out)   :: conv_flag !< Convergence flag: 1=converged, 0=not
  integer,                intent(out)   :: iters !< The number of iterations used
  integer,                intent(in)    :: Is_sum !< Starting i-index for global sums
  integer,                intent(in)    :: Js_sum !< Starting j-index for global sums
  integer,                intent(in)    :: Ie_sum !< Ending i-index for global sums
  integer,                intent(in)    :: Je_sum !< Ending j-index for global sums
  integer,                intent(in)    :: Iscq_sv !< Starting i-index for sum_vec arrays
  integer,                intent(in)    :: Jscq_sv !< Starting j-index for sum_vec arrays
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          optional, intent(in) :: h_shelf !< Ice shelf thickness on tracer grid [Z ~> m]

  real, dimension(SZDIB_(G),SZDJB_(G))  :: u_curr  !< Frozen current iterate u^k, used to evaluate basal friction
                                                   !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G))  :: v_curr  !< Frozen current iterate v^k, used to evaluate basal friction
                                                   !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) ::  &
                        Ru, Rv, &         ! Residuals (r) [R L3 Z T-2 ~> m kg s-2]
                        Zu, Zv, &         ! Preconditioned residuals (z = M^-1 r) [L T-1 ~> m s-1]
                        Du, Dv, &         ! Search directions (p) [L T-1 ~> m s-1]
                        Qu, Qv            ! A * p [R L3 Z T-2 ~> m kg s-2]
  real, dimension(SZDIB_(G),SZDJB_(G),2) :: sum_vec_3d ! Pointwise products for global sums.
                                ! sum_vec_3d(:,:,1): r^2 [kg2 m2 s-4] or z·q [kg m2 s-3] (context-dependent)
                                ! sum_vec_3d(:,:,2): z·w or q·(M^-1 q) [kg m2 s-3]
  real    :: alpha        ! Step length [nondim]
  real    :: beta         ! Direction update coefficient [nondim]
  real    :: r_norm_sq    ! Squared residual norm [kg2 m2 s-4]
  real    :: z_w_sum      ! Inner product (z_k, A z_k); beta denominator [kg m2 s-3]
  real    :: z_w_sum_new  ! Inner product (z_{k+1}, A z_{k+1}); beta numerator [kg m2 s-3]
  real    :: z_q_sum      ! Inner product (z_k, A p_k); alpha numerator [kg m2 s-3]
  real    :: q_s_sum      ! Inner product (A p_k, M^-1 A p_k); alpha denom [kg m2 s-3]
  real    :: resid0tol2   ! Convergence threshold: tol^2 * ||r_0||^2 [kg2 m2 s-4]
  real    :: sv3dsum      ! Unused scalar return from reproducing_sum [various]
  real    :: sv3dsums(2)  ! Component sums from reproducing_sum
                          ! sv3dsums(1): r^2 or z·q [kg2 m2 s-4 or kg m2 s-3] (context-dependent)
                          ! sv3dsums(2): z·w or q·M^-1 q [kg m2 s-3]
  real    :: resid2_scale ! Scaling for squared-stress inner products [T4 kg2 m2 R-2 Z-2 L-6 s-4 ~> 1]
  integer :: iter, i, j, isc, iec, jsc, jec
  integer :: Isdq, Iedq, Jsdq, Jedq, Iscq, Iecq, Jscq, Jecq

  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  resid2_scale = ((US%RZ_to_kg_m2*US%L_to_m)*US%L_T_to_m_s**2)**2

  ! Initialize CR-specific arrays
  Ru(:,:) = 0 ; Rv(:,:) = 0 ; Zu(:,:) = 0 ; Zv(:,:) = 0
  Du(:,:) = 0 ; Dv(:,:) = 0 ; Qu(:,:) = 0 ; Qv(:,:) = 0

  ! r_0 = b - A*x_0
  Ru(:,:) = (RHSu(:,:) - Au(:,:)) ; Rv(:,:) = (RHSv(:,:) - Av(:,:))

  ! current velocities used in CG_action for basal drag
  u_curr(:,:) = u_shlf(:,:) ; v_curr(:,:) = v_shlf(:,:)

  ! z_0 = M^-1 r_0
  do J=Jsdq,Jedq ; do I=Isdq,Iedq
     if (CS%umask(I,J) == 1) Zu(I,J) = Ru(I,J) * IDIAGu(I,J)
     if (CS%vmask(I,J) == 1) Zv(I,J) = Rv(I,J) * IDIAGv(I,J)
  enddo ; enddo

  ! p_0 = z_0
  Du(:,:) = Zu(:,:) ; Dv(:,:) = Zv(:,:)

  ! Compute A * z_0
  Au(:,:) = 0 ; Av(:,:) = 0
  call CG_action(CS, Au, Av, Zu, Zv, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                 H_node, CS%ice_visc, CS%bed_elev, u_curr, v_curr, &
                 G, US, isc-1, iec+1, jsc-1, jec+1, rhoi_rhow, h_shelf=h_shelf)
  call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE)

  ! q_0 = A * p_0
  Qu(:,:) = Au(:,:) ; Qv(:,:) = Av(:,:)

  ! Initial Norms
  sum_vec_3d(:,:,:) = 0.0 ; sv3dsums(1:2) = 0.0
  do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
    if (CS%umask(I,J) == 1) then
      sum_vec_3d(I,J,1) = resid2_scale * Ru(I,J)**2
      sum_vec_3d(I,J,2) = resid_scale  * (Zu(I,J) * Au(I,J))
    endif
    if (CS%vmask(I,J) == 1) then
      sum_vec_3d(I,J,1) = sum_vec_3d(I,J,1) + resid2_scale * Rv(I,J)**2
      sum_vec_3d(I,J,2) = sum_vec_3d(I,J,2) + resid_scale  * (Zv(I,J) * Av(I,J))
    endif
  enddo ; enddo
  sv3dsum = reproducing_sum( sum_vec_3d(:,:,1:2), Is_sum, Ie_sum, Js_sum, Je_sum, sums=sv3dsums(1:2) )

  r_norm_sq = sv3dsums(1)
  z_w_sum   = sv3dsums(2)

  resid0tol2 = CS%cg_tol_current**2 * r_norm_sq
  conv_flag = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                                !!
  !! MAIN CONJUGATE RESIDUAL LOOP   !!
  !!                                !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do iter = 1, CS%cg_max_iterations

    ! --- STEP 1: alpha = (z_k, q_k) / (q_k, M^-1 q_k) ---
    sum_vec_3d(:,:,:) = 0.0 ; sv3dsums(1:2) = 0.0
    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
       if (CS%umask(I,J) == 1) then
           sum_vec_3d(I,J,1) = resid_scale * (Zu(I,J) * Qu(I,J))
           ! Order matters to prevent float overflow: Q * (Q * IDiag)
           sum_vec_3d(I,J,2) = resid_scale * (Qu(I,J) * (Qu(I,J) * IDIAGu(I,J)))
       endif
       if (CS%vmask(I,J) == 1) then
           sum_vec_3d(I,J,1) = sum_vec_3d(I,J,1) + resid_scale * (Zv(I,J) * Qv(I,J))
           sum_vec_3d(I,J,2) = sum_vec_3d(I,J,2) + resid_scale * (Qv(I,J) * (Qv(I,J) * IDIAGv(I,J)))
       endif
    enddo ; enddo
    sv3dsum = reproducing_sum( sum_vec_3d(:,:,1:2), Is_sum, Ie_sum, Js_sum, Je_sum, sums=sv3dsums(1:2) )

    z_q_sum = sv3dsums(1)
    q_s_sum = sv3dsums(2)

    if (q_s_sum == 0.0) then
      iters = iter
      conv_flag = 1
      exit
    endif
    alpha = z_q_sum / q_s_sum

    ! --- STEP 2: Update x, r, and z (Fused over Full Domain) ---
    ! Zu halos are populated here since the loop covers Jsdq..Jedq; no pass_vector needed.
    do J=Jsdq,Jedq ; do I=Isdq,Iedq
       if (CS%umask(I,J) == 1) then
           u_shlf(I,J) = u_shlf(I,J) + alpha * Du(I,J)
           Ru(I,J) = Ru(I,J) - alpha * Qu(I,J)
           Zu(I,J) = Ru(I,J) * IDIAGu(I,J)
       endif
       if (CS%vmask(I,J) == 1) then
           v_shlf(I,J) = v_shlf(I,J) + alpha * Dv(I,J)
           Rv(I,J) = Rv(I,J) - alpha * Qv(I,J)
           Zv(I,J) = Rv(I,J) * IDIAGv(I,J)
       endif
    enddo ; enddo

    ! --- STEP 3: w_{k+1} = A z_{k+1} ---
    Au(:,:) = 0 ; Av(:,:) = 0
    call CG_action(CS, Au, Av, Zu, Zv, Phi, Phisub, CS%umask, CS%vmask, hmask, &
                   H_node, CS%ice_visc, CS%bed_elev, u_curr, v_curr, &
                   G, US, isc-1, iec+1, jsc-1, jec+1, rhoi_rhow, h_shelf=h_shelf)
    call pass_vector(Au, Av, G%domain, TO_ALL, BGRID_NE)

    ! --- STEP 4: beta and convergence check ---
    sum_vec_3d(:,:,:) = 0.0 ; sv3dsums(1:2) = 0.0
    do J=Jscq_sv,Jecq ; do I=Iscq_sv,Iecq
       if (CS%umask(I,J) == 1) then
           sum_vec_3d(I,J,1) = resid2_scale * Ru(I,J)**2
           sum_vec_3d(I,J,2) = resid_scale  * (Zu(I,J) * Au(I,J))
       endif
       if (CS%vmask(I,J) == 1) then
           sum_vec_3d(I,J,1) = sum_vec_3d(I,J,1) + resid2_scale * Rv(I,J)**2
           sum_vec_3d(I,J,2) = sum_vec_3d(I,J,2) + resid_scale  * (Zv(I,J) * Av(I,J))
       endif
    enddo ; enddo
    sv3dsum = reproducing_sum( sum_vec_3d(:,:,1:2), Is_sum, Ie_sum, Js_sum, Je_sum, sums=sv3dsums(1:2) )

    r_norm_sq = sv3dsums(1)
    z_w_sum_new = sv3dsums(2)

    if (r_norm_sq <= resid0tol2 .or. z_w_sum==0.0) then
      iters = iter
      conv_flag = 1
      exit
    endif

    beta = z_w_sum_new / z_w_sum
    z_w_sum = z_w_sum_new

    ! --- STEP 5: Update p and q ---
    do J=Jsdq,Jedq ; do I=Isdq,Iedq
       if (CS%umask(I,J) == 1) then
           Du(I,J) = Zu(I,J) + beta * Du(I,J)
           Qu(I,J) = Au(I,J) + beta * Qu(I,J)
       endif
       if (CS%vmask(I,J) == 1) then
           Dv(I,J) = Zv(I,J) + beta * Dv(I,J)
           Qv(I,J) = Av(I,J) + beta * Qv(I,J)
       endif
    enddo ; enddo

  enddo ! end of CR loop

end subroutine ice_shelf_solve_inner_CR

subroutine ice_shelf_advect_thickness_x(CS, G, LB, time_step, hmask, h0, h_after_uflux, uh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(loop_bounds_type), intent(in)    :: LB   !< Loop bounds structure.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_uflux !< The ice shelf thicknesses after
                                              !! the zonal mass fluxes [Z ~> m].
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization


  integer :: i, j
  integer :: ish, ieh, jsh, jeh
  real :: u_face     ! Zonal velocity at a face [L T-1 ~> m s-1]
  real :: h_face     ! Thickness at a face for transport [Z ~> m]
  real :: slope_lim  ! The value of the slope limiter, in the range of 0 to 2 [nondim]

!  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec
!  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  ! hmask coded values: 1) fully covered; 2) partly covered - no export; 3) Specified boundary condition
  ! relevant u_face_mask coded values: 1) Normal interior point; 4) Specified flux BC

  do j=jsh,jeh ; do I=ish-1,ieh
    if (CS%u_face_mask(I,j) == 4.) then ! The flux itself is a specified boundary condition.
      uh_ice(I,j) = (time_step * G%dyCu(I,j)) * CS%u_flux_bdry_val(I,j)
    elseif ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .or. (hmask(i+1,j) == 1 .or. hmask(i+1,j) == 3)) then
      u_face = 0.5 * (CS%u_shelf(I,J-1) + CS%u_shelf(I,J))
      h_face = 0.0 ! This will apply when the source cell is iceless or not fully ice covered.

      if (u_face > 0) then
        if (hmask(i,j) == 3) then ! This is a open boundary inflow from the west
          h_face = CS%h_bdry_val(i,j)
        elseif (hmask(i,j) == 1) then ! There can be eastward flow through this face.
          if ((hmask(i-1,j) == 1 .or. hmask(i-1,j) == 3) .and. &
            (hmask(i+1,j) == 1 .or. hmask(i+1,j) == 3)) then
            slope_lim = slope_limiter(h0(i,j)-h0(i-1,j), h0(i+1,j)-h0(i,j))
            ! This is a 2nd-order centered scheme with a slope limiter.  We could try PPM here.
            h_face = h0(i,j) - slope_lim * (0.5 * (h0(i,j)-h0(i+1,j)))
            if (associated(CS%phi_x_FV)) CS%phi_x_FV(I,j) = slope_lim
          else
            h_face = h0(i,j)
          endif
        endif
      else
        if (hmask(i+1,j) == 3) then ! This is a open boundary inflow from the east
          h_face = CS%h_bdry_val(i+1,j)
        elseif (hmask(i+1,j) == 1) then
          if ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .and. &
            (hmask(i+2,j) == 1 .or. hmask(i+2,j) == 3)) then
            slope_lim = slope_limiter(h0(i+1,j)-h0(i,j), h0(i+2,j)-h0(i+1,j))
            h_face = h0(i+1,j) - slope_lim * (0.5 * (h0(i+2,j)-h0(i+1,j)))
            if (associated(CS%phi_x_FV)) CS%phi_x_FV(I,j) = slope_lim
          else
            h_face = h0(i+1,j)
          endif
        endif
      endif

      uh_ice(I,j) = (time_step * G%dyCu(I,j)) * (u_face * h_face)
    else
      uh_ice(I,j) = 0.0
    endif
  enddo ; enddo

  do j=jsh,jeh ; do i=ish,ieh
    if (hmask(i,j) /= 3) &
      h_after_uflux(i,j) = h0(i,j) + (uh_ice(I-1,j) - uh_ice(I,j)) * G%IareaT(i,j)

     ! Update the masks of cells that have gone from no ice to partial ice.
    if ((hmask(i,j) == 0) .and. ((uh_ice(I-1,j) > 0.0) .or. (uh_ice(I,j) < 0.0))) hmask(i,j) = 2
  enddo ; enddo

end subroutine ice_shelf_advect_thickness_x

subroutine ice_shelf_advect_thickness_y(CS, G, LB, time_step, hmask, h0, h_after_vflux, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(loop_bounds_type), intent(in)    :: LB !< Loop bounds structure.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                              !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_vflux !< The ice shelf thicknesses after
                                              !! the meridional mass fluxes [Z ~> m].
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization


  integer :: i, j
  integer :: ish, ieh, jsh, jeh
  real :: v_face     ! Pseudo-meridional velocity at a face [L T-1 ~> m s-1]
  real :: h_face     ! Thickness at a face for transport [Z ~> m]
  real :: slope_lim  ! The value of the slope limiter, in the range of 0 to 2 [nondim]

  ish = LB%ish ; ieh = LB%ieh ; jsh = LB%jsh ; jeh = LB%jeh

  ! hmask coded values: 1) fully covered; 2) partly covered - no export; 3) Specified boundary condition
  ! relevant u_face_mask coded values: 1) Normal interior point; 4) Specified flux BC

  do J=jsh-1,jeh ; do i=ish,ieh
    if (CS%v_face_mask(i,J) == 4.) then ! The flux itself is a specified boundary condition.
      vh_ice(i,J) = (time_step * G%dxCv(i,J)) * CS%v_flux_bdry_val(i,J)
    elseif ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .or. (hmask(i,j+1) == 1 .or. hmask(i,j+1) == 3)) then
      v_face = 0.5 * (CS%v_shelf(I-1,J) + CS%v_shelf(I,J))
      h_face = 0.0 ! This will apply when the source cell is iceless or not fully ice covered.

      if (v_face > 0) then
        if (hmask(i,j) == 3) then ! This is a open boundary inflow from the south
          h_face = CS%h_bdry_val(i,j)
        elseif (hmask(i,j) == 1) then ! There can be northward flow through this face.
          if ((hmask(i,j-1) == 1 .or. hmask(i,j-1) == 3) .and. &
            (hmask(i,j+1) == 1 .or. hmask(i,j+1) == 3)) then
            slope_lim = slope_limiter(h0(i,j)-h0(i,j-1), h0(i,j+1)-h0(i,j))
            ! This is a 2nd-order centered scheme with a slope limiter.  We could try PPM here.
            h_face = h0(i,j) - slope_lim * (0.5 * (h0(i,j)-h0(i,j+1)))
            if (associated(CS%phi_y_FV)) CS%phi_y_FV(i,J) = slope_lim
          else
            h_face = h0(i,j)
          endif
        endif
      else
        if (hmask(i,j+1) == 3) then ! This is a open boundary inflow from the north
          h_face = CS%h_bdry_val(i,j+1)
        elseif (hmask(i,j+1) == 1) then
          if ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .and. &
            (hmask(i,j+2) == 1 .or. hmask(i,j+2) == 3)) then
            slope_lim = slope_limiter(h0(i,j+1)-h0(i,j), h0(i,j+2)-h0(i,j+1))
            h_face = h0(i,j+1) - slope_lim * (0.5 * (h0(i,j+2)-h0(i,j+1)))
            if (associated(CS%phi_y_FV)) CS%phi_y_FV(i,J) = slope_lim
          else
            h_face = h0(i,j+1)
          endif
        endif
      endif

      vh_ice(i,J) = (time_step * G%dxCv(i,J)) * (v_face * h_face)
    else
      vh_ice(i,J) = 0.0
    endif
  enddo ; enddo

  do j=jsh,jeh ; do i=ish,ieh
    if (hmask(i,j) /= 3) &
      h_after_vflux(i,j) = h0(i,j) + (vh_ice(i,J-1) - vh_ice(i,J)) * G%IareaT(i,j)

    ! Update the masks of cells that have gone from no ice to partial ice.
    if ((hmask(i,j) == 0) .and. ((vh_ice(i,J-1) > 0.0) .or. (vh_ice(i,J) < 0.0))) hmask(i,j) = 2
  enddo ; enddo

end subroutine ice_shelf_advect_thickness_y

subroutine shelf_advance_front(CS, ISS, G, hmask, uh_ice, vh_ice, calving)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                           !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  logical, optional, intent(in)         :: calving !< If True, shelf_advance front is being
                                                                 !! used for calving from a static ice front
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                              !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]

  ! in this subroutine we go through the computational cells only and, if they are empty or partial cells,
  ! we find the reference thickness and update the shelf mass and partial area fraction and the hmask if necessary

  ! if any cells go from partial to complete, we then must set the thickness, update hmask accordingly,
  ! and divide the overflow across the adjacent EMPTY (not partly-covered) cells.
  ! (it is highly unlikely there will not be any; in which case this will need to be rethought.)

  ! most likely there will only be one "overflow". If not, though, a pass_var of all relevant variables
  ! is done; there will therefore be a loop which, in practice, will hopefully not have to go through
  ! many iterations

  ! when 3d advected scalars are introduced, they will be impacted by what is done here

  ! flux_enter(isd:ied,jsd:jed,1:4): if cell is not ice-covered, gives flux of ice into cell from kth boundary
  !
  !   from eastern neighbor:  flux_enter(:,:,1)
  !   from western neighbor:  flux_enter(:,:,2)
  !   from southern neighbor: flux_enter(:,:,3)
  !   from northern neighbor: flux_enter(:,:,4)
  !
  !        o--- (4) ---o
  !        |           |
  !       (1)         (2)
  !        |           |
  !        o--- (3) ---o
  !

  integer :: i, j, isc, iec, jsc, jec, n_flux, k, iter_count
  integer :: i_off, j_off
  integer :: iter_flag

  real :: h_reference ! A reference thicknesss based on neighboring cells [Z ~> m]
  real :: h_reference_ew !contribution to reference thickness from east + west cells [Z ~> m]
  real :: h_reference_ns !contribution to reference thickness from north + south cells [Z ~> m]
  real :: tot_flux    ! The total ice mass flux [Z L2 ~> m3]
  real :: tot_flux_ew ! The contribution to total ice mass flux from east + west cells [Z L2 ~> m3]
  real :: tot_flux_ns ! The contribution to total ice mass flux from north + south cells [Z L2 ~> m3]
  real :: partial_vol ! The volume covered by ice shelf [Z L2 ~> m3]
  real :: dxdyh       ! Cell area [L2 ~> m2]
  logical :: ice_shelf_calving ! True if this subroutine is being used for ice-shelf calving
  character(len=160) :: mesg  ! The text of an error message
  integer, dimension(4) :: mapi, mapj, new_partial
  real, dimension(SZDI_(G),SZDJ_(G),4) :: flux_enter  ! The ice volume flux into the
                                              ! cell through the 4 cell boundaries [Z L2 ~> m3].
  real, dimension(SZDI_(G),SZDJ_(G),4) :: flux_enter_replace ! An updated ice volume flux into the
                                              ! cell through the 4 cell boundaries [Z L2 ~> m3].

  if (present(calving)) then
    ice_shelf_calving = calving
  else
    ice_shelf_calving = .false.
  endif

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = G%idg_offset ; j_off = G%jdg_offset
  iter_count = 0 ; iter_flag = 1

  flux_enter(:,:,:) = 0.0
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if ((hmask(i,j) == 0) .or. (hmask(i,j) == 2)) then
      flux_enter(i,j,1) = max(uh_ice(I-1,j), 0.0)
      flux_enter(i,j,2) = max(-uh_ice(I,j), 0.0)
      flux_enter(i,j,3) = max(vh_ice(i,J-1), 0.0)
      flux_enter(i,j,4) = max(-vh_ice(i,J), 0.0)
    endif
  enddo ; enddo

  mapi(1) = -1 ; mapi(2) = 1 ; mapi(3:4) = 0
  mapj(3) = -1 ; mapj(4) = 1 ; mapj(1:2) = 0

  do while (iter_flag == 1)

    iter_flag = 0

    if (iter_count > 0) then
      flux_enter(:,:,:) = flux_enter_replace(:,:,:)
    endif
    flux_enter_replace(:,:,:) = 0.0

    iter_count = iter_count + 1

    ! if iter_count >= 3 then some halo updates need to be done...
    if (iter_count==3) then
      call MOM_error(FATAL, "MOM_ice_shelf_dyn.F90, shelf_advance_front iter >=3.")
    endif

    do j=jsc-1,jec+1

      if (CS%reentrant_y .OR. (((j+j_off) <= G%domain%njglobal) .AND. &
          ((j+j_off) >= 1))) then

        do i=isc-1,iec+1

          if (CS%reentrant_x .OR. (((i+i_off) <= G%domain%niglobal) .AND. &
              ((i+i_off) >= 1))) then
            ! first get reference thickness by averaging over cells that are fluxing into this cell
            n_flux = 0
            h_reference_ew = 0.0
            h_reference_ns = 0.0
            tot_flux_ew = 0.0
            tot_flux_ns = 0.0

            do k=1,2
              if (flux_enter(i,j,k) > 0) then
                n_flux = n_flux + 1
                h_reference_ew = h_reference_ew + flux_enter(i,j,k) * ISS%h_shelf(i+2*k-3,j)
                !h_reference = h_reference + ISS%h_shelf(i+2*k-3,j)
                tot_flux_ew = tot_flux_ew + flux_enter(i,j,k)
                flux_enter(i,j,k) = 0.0
              endif
            enddo

            do k=1,2
              if (flux_enter(i,j,k+2) > 0) then
                n_flux = n_flux + 1
                h_reference_ns = h_reference_ns + flux_enter(i,j,k+2) * ISS%h_shelf(i,j+2*k-3)
                !h_reference = h_reference + ISS%h_shelf(i,j+2*k-3)
                tot_flux_ns = tot_flux_ns + flux_enter(i,j,k+2)
                flux_enter(i,j,k+2) = 0.0
              endif
            enddo

            h_reference = h_reference_ew + h_reference_ns
            tot_flux = tot_flux_ew + tot_flux_ns

            if (n_flux > 0) then
              dxdyh = G%areaT(i,j)
              h_reference = h_reference / tot_flux
              !h_reference = h_reference / real(n_flux)
              partial_vol = ISS%h_shelf(i,j) * ISS%area_shelf_h(i,j) + tot_flux
              ! The partial-fill overwrites h_shelf with the donor cell mean;
              ! set the DG nodal field to that same donor mean so
              ! nodal_cell_mean matches ISS%h_shelf. Leaving the corners at
              ! 0 (the prior code path) collapses Hmin_B at this cell's 4
              ! B-nodes to 0, which lets neighbour ice cells show corner
              ! jumps as large as their own Hbar through the nodal limiter
              ! envelope.
              if (CS%use_DG_thickness) then
                CS%h_nodal(i,j,:,:) = h_reference
              endif

              if (ice_shelf_calving) then
                !Mark calving cells as "underfilled" (hmask = 2), even in the unlikely case that they
                !are exactly covered or overflowed. Hmask, area_shelf_h and h_shelf will be reset to zero
                !after being used to calculate the calving mass for the cell.
                ISS%hmask(i,j) = 2
                ISS%area_shelf_h(i,j) = partial_vol / h_reference
                ISS%h_shelf(i,j) = h_reference
              elseif ((partial_vol / G%areaT(i,j)) == h_reference) then ! cell is exactly covered, no overflow
                if (ISS%hmask(i,j)/=3) ISS%hmask(i,j) = 1
                ISS%h_shelf(i,j) = h_reference
                ISS%area_shelf_h(i,j) = G%areaT(i,j)
              elseif ((partial_vol / G%areaT(i,j)) < h_reference) then
                ISS%hmask(i,j) = 2
               !  ISS%mass_shelf(i,j) = partial_vol * CS%density_ice
                ISS%area_shelf_h(i,j) = partial_vol / h_reference
                ISS%h_shelf(i,j) = h_reference
              else

                if (ISS%hmask(i,j)/=3) ISS%hmask(i,j) = 1
                ISS%area_shelf_h(i,j) = G%areaT(i,j)
                !h_temp(i,j) = h_reference
                partial_vol = partial_vol - h_reference * G%areaT(i,j)

                iter_flag  = 1

                n_flux = 0 ; new_partial(:) = 0

                do k=1,2
                  if (CS%u_face_mask(I-2+k,j) == 2) then
                    n_flux = n_flux + 1
                  elseif (ISS%hmask(i+2*k-3,j) == 0) then
                    n_flux = n_flux + 1
                    new_partial(k) = 1
                  endif
                  if (CS%v_face_mask(i,J-2+k) == 2) then
                    n_flux = n_flux + 1
                  elseif (ISS%hmask(i,j+2*k-3) == 0) then
                    n_flux = n_flux + 1
                    new_partial(k+2) = 1
                  endif
                enddo

                if (n_flux == 0) then ! there is nowhere to put the extra ice!
                  ISS%h_shelf(i,j) = h_reference + partial_vol / G%areaT(i,j)
                else
                  ISS%h_shelf(i,j) = h_reference

                  do k=1,2
                    if (new_partial(k) == 1) &
                      flux_enter_replace(i+2*k-3,j,3-k) = partial_vol / real(n_flux)
                    if (new_partial(k+2) == 1) &
                      flux_enter_replace(i,j+2*k-3,5-k) = partial_vol / real(n_flux)
                  enddo
                endif

              endif ! Parital_vol test.
            endif ! n_flux gt 0 test.

          endif
        enddo ! j-loop
      endif
    enddo

  !  call max_across_PEs(iter_flag)

  enddo ! End of do while(iter_flag) loop

  call max_across_PEs(iter_count)

  if (is_root_pe() .and. (iter_count > 1)) then
    write(mesg,*) "shelf_advance_front: ", iter_count, " max iterations"
    call MOM_mesg(mesg, 5)
  endif

end subroutine shelf_advance_front

!> Calculate total horizontal flux in/out of the domain. This subroutine could be used to  calculate the
!! stocks in the hole that may appear grid at the South Pole. The flux is calculated over edges of the
!! computational domain with non-zero boundary conditions set for velocity or flux.
subroutine calculate_flux_inout(CS, ISS, G, uh_ice, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(inout) :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJ_(G)), intent(in) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G)), intent(in) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]
  integer :: i, j, isc, iec, jsc, jec
  integer :: i_off, j_off
  real, dimension(SZDIB_(G),SZDJ_(G)) :: u_flux ! Accumulated zonal flux in/out of the domain
                                                ! (outward is positive) [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G)) :: v_flux ! Accumulated meridional flux in/out of the domain
                                                  ! (outward is positive) [Z L2 ~> m3]
  integer :: Isdq, Iedq, Jsdq, Jedq
  integer :: Iscq, Iecq, Jscq, Jecq
  integer :: Is_sum, Js_sum, Ie_sum, Je_sum ! Loop bounds for global sums or arrays starting at 1.
  integer :: Iscq_sv, Jscq_sv ! Starting loop bound for sum_vec
  character(len=160) :: mesg  ! The text of an error message

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  Isdq = G%IsdB ; Iedq = G%IedB ; Jsdq = G%JsdB ; Jedq = G%JedB
  Iscq = G%IscB ; Iecq = G%IecB ; Jscq = G%JscB ; Jecq = G%JecB
  i_off = G%idg_offset ; j_off = G%jdg_offset

  u_flux(:,:) = 0.0
  v_flux(:,:) = 0.0

  !Southern boundary (symmetric only, flux here is zero if non-symmetric)
  if (G%symmetric .and. (.not. CS%reentrant_y)) then
    if (jsc+j_off == G%jsg) then
      J=jsc-1
      do i = isc,iec
        if (CS%v_face_mask(i,J)==3 .or. CS%v_face_mask(i,J)==4 .or. CS%v_face_mask(i,J)==6) then
          v_flux(i,J) = -vh_ice(i,J)
        endif
      enddo
    endif
  endif

  !Western boundary (symmetric only, flux here is zero if non-symmetric)
  if (G%symmetric .and. (.not. CS%reentrant_x)) then
    if (isc+i_off == G%isg) then
      I=isc-1
      do j = jsc,jec
        if (CS%u_face_mask(I,j)>=3 .and. CS%u_face_mask(I,j)<=5) then
          u_flux(I,j) = -uh_ice(I,j)
        endif
      enddo
    endif
  endif

  !Northern boundary
  if (.not. CS%reentrant_y) then
    if (jec+j_off == G%domain%njglobal) then
      J=jec
      do i = isc,iec
        if (CS%v_face_mask(i,J)==3 .or. CS%v_face_mask(i,J)==4 .or. CS%v_face_mask(i,J)==6) then
          v_flux(i,J) = vh_ice(i,J)
        endif
      enddo
    endif
  endif

  !Eastern boundary
  if (.not. CS%reentrant_x) then
    if (iec+i_off == G%domain%niglobal) then
      I=iec
      do j = jsc,jec
        if (CS%u_face_mask(I,j)>=3 .and. CS%u_face_mask(I,j)<=5) then
          u_flux(I,j) = uh_ice(I,j)
        endif
      enddo
    endif
  endif

  ! Determine the loop limits for sums, bearing in mind that the arrays will be starting at 1.
  ! Includes the edge of the tile is at the western/southern bdry (if symmetric)
  if ((isc+G%idg_offset==G%isg) .and. (.not. CS%reentrant_x)) then
    Is_sum = Iscq + (1-Isdq) ; Iscq_sv = Iscq
  else
    Is_sum = isc  + (1-Isdq) ; Iscq_sv = isc
  endif
  if ((jsc+G%jdg_offset==G%jsg) .and. (.not. CS%reentrant_y)) then
    Js_sum = Jscq + (1-Jsdq) ; Jscq_sv = Jscq
  else
    Js_sum = jsc + (1-Jsdq) ; Jscq_sv = jsc
  endif
  Ie_sum = Iecq + (1-Isdq) ; Je_sum = Jecq + (1-Jsdq)

  !Total accumulated flux in/out of the domain edges (outward is positive)
  ISS%tot_flux_inout = reproducing_sum(u_flux, Is_sum, Ie_sum, Js_sum, Je_sum) + &
                       reproducing_sum(v_flux, Is_sum, Ie_sum, Js_sum, Je_sum)

  write(mesg,*) 'Flux in/out of domain', ISS%tot_flux_inout * CS%density_ice * G%US%RZL2_to_kg
  call MOM_mesg("MOM6-IS: "//trim(mesg))
end subroutine calculate_flux_inout

!> Apply a very simple calving law using a minimum thickness rule
subroutine ice_shelf_min_thickness_calve(G, h_shelf, area_shelf_h, hmask, thickness_calve, halo, h_nodal)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                  intent(in)    :: thickness_calve !< The thickness at which to trigger calving [Z ~> m].
  integer,     optional, intent(in)    :: halo  !< The number of halo points to use.  If not present,
                                                !! work on the entire data domain.
  real, dimension(SZDI_(G),SZDJ_(G),2,2), optional, intent(inout) :: h_nodal !< Nodal Q1 thickness [Z ~> m],
                                                                    !! zeroed when a cell calves.
  integer :: i, j, is, ie, js, je

  if (present(halo)) then
    is = G%isc - halo ; ie = G%iec + halo ; js = G%jsc - halo ; je = G%jec + halo
  else
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif

  do j=js,je ; do i=is,ie
    if ((h_shelf(i,j) < thickness_calve) .and. (area_shelf_h(i,j) > 0.)) then
      h_shelf(i,j) = 0.0
      area_shelf_h(i,j) = 0.0
      hmask(i,j) = 0.0
      if (present(h_nodal)) h_nodal(i,j,:,:) = 0.0
    endif
  enddo ; enddo

end subroutine ice_shelf_min_thickness_calve

subroutine calve_to_mask(G, h_shelf, area_shelf_h, hmask, calve_mask, h_nodal)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in)    :: calve_mask !< A mask that indicates where the ice
                                                             !! shelf can exist, and where it will calve.
  real, dimension(SZDI_(G),SZDJ_(G),2,2), optional, intent(inout) :: h_nodal !< Nodal Q1 thickness [Z ~> m],
                                                                    !! zeroed when a cell calves.

  integer                        :: i,j

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((calve_mask(i,j) == 0.0) .and. (hmask(i,j) /= 0.0)) then
      h_shelf(i,j) = 0.0
      area_shelf_h(i,j) = 0.0
      hmask(i,j) = 0.0
      if (present(h_nodal)) h_nodal(i,j,:,:) = 0.0
    endif
  enddo ; enddo

end subroutine calve_to_mask

!> Calculate driving stress using cell-centered bed elevation and ice thickness
subroutine calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, OD)
  type(ice_shelf_dyn_CS), intent(in)   :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(in)    :: ISS !< A structure with elements that describe
                                             !! the ice-shelf state
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: OD  !< ocean floor depth at tracer points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudx  !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudy  !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]


! driving stress!

! ! taudx and taudy will hold driving stress in the x- and y- directions when done.
!    they will sit on the BGrid, and so their size depends on whether the grid is symmetric
!
! Since this is a finite element solve, they will actually have the form \int \Phi_i rho g h \nabla s
!
! OD -this is important and we do not yet know where (in MOM) it will come from. It represents
!     "average" ocean depth -- and is needed to find surface elevation
!    (it is assumed that base_ice = bed + OD)

  real, dimension(SIZE(OD,1),SIZE(OD,2))  :: S     ! surface elevation [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)) :: sx_e, sy_e !element contributions to driving stress
  real    :: rho, rhow, rhoi_rhow ! Ice and ocean densities [R ~> kg m-3]
  real    :: sx, sy    ! Ice shelf top slopes at tracer points [Z L-1 ~> nondim]
  real    :: neumann_val ! [R Z L2 T-2 ~> kg s-2]
  real    :: grav      ! The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real    :: scale     ! Scaling factor used to ensure surface slope magnitude does not exceed CS%max_surface_slope
  logical :: valid_N, valid_S, valid_E, valid_W
  integer :: i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, is, js, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, isc, jsc, iec, jec
  integer :: i_off, j_off

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
!  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
!  iegq = G%iegB ; jegq = G%jegB
!  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  gisc = 1 ; gjsc = 1
!  giec = G%domain%niglobal+G%domain%nihalo ; gjec = G%domain%njglobal+G%domain%njhalo
  giec = G%domain%niglobal ; gjec = G%domain%njglobal
!  is = iscq - 1 ; js = jscq - 1
  i_off = G%idg_offset ; j_off = G%jdg_offset


  rho =  CS%density_ice
  rhow = CS%density_ocean_avg
  grav = CS%g_Earth
  rhoi_rhow = rho/rhow
  ! prelim - go through and calculate S

  if (CS%GL_couple) then
    do j=jsc-2,jec+2 ; do i=isc-2,iec+2
      S(i,j) = -CS%bed_elev(i,j) + (OD(i,j) + max(ISS%h_shelf(i,j),CS%min_h_shelf))
    enddo ; enddo
  else
    ! check whether the ice is floating or grounded
    do j=jsc-2,jec+2 ; do i=isc-2,iec+2
      if (rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf) - CS%bed_elev(i,j) <= 0) then
        S(i,j) = (1 - rhoi_rhow)*max(ISS%h_shelf(i,j),CS%min_h_shelf)
      else
        S(i,j) = max(ISS%h_shelf(i,j),CS%min_h_shelf)-CS%bed_elev(i,j)
      endif
    enddo ; enddo
  endif

  call pass_var(S, G%domain)

  do j=jsc-1,jec+1
    do i=isc-1,iec+1

      if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then
        ! we are inside the global computational bdry, at an ice-filled cell

        ! Calculate the x-direction surface slope at tracer points.
        sx = 0.0
        valid_E = (ISS%hmask(i+1,j) == 1 .or. ISS%hmask(i+1,j) == 3)
        valid_W = (ISS%hmask(i-1,j) == 1 .or. ISS%hmask(i-1,j) == 3)
        if (CS%shelf_top_slope_bugs) then
          if (((i+i_off) == gisc) .and. (.not.CS%reentrant_x)) then ! at west computational bdry
            if (valid_E) sx = (S(i+1,j)-S(i,j)) / G%dxT(i,j)
          elseif (((i+i_off) == giec) .and. (.not.CS%reentrant_x)) then ! at east computational bdry
            if (valid_W) sx = (S(i,j)-S(i-1,j)) / G%dxT(i,j)
          elseif (valid_E .and. valid_W) then
            ! This is the usual interior point
            sx = (S(i+1,j) - S(i-1,j)) / (G%dxT(i,j) + G%dxT(i-1,j))
          elseif (valid_E) then
            sx = (S(i+1,j) - S(i,j)) / (G%dxT(i,j) + G%dxT(i+1,j))
          elseif (valid_W) then
            sx = (S(i,j) - S(i-1,j)) / (G%dxT(i,j) + G%dxT(i-1,j))
          endif
        else ! Correct the bugs in the version above.
          if (((i+i_off) == gisc) .and. (.not.CS%reentrant_x)) then ! at west computational bdry
            if (valid_E) sx = (S(i+1,j) - S(i,j)) * G%IdxCu(I,j)
          elseif (((i+i_off) == giec) .and. (.not.CS%reentrant_x)) then ! at east computational bdry
            if (valid_W) sx = (S(i,j) - S(i-1,j)) * G%IdxCu(I-1,j)
          elseif (valid_E .and. valid_W) then
            ! This is the usual interior point
            sx = 0.5*(S(i+1,j) - S(i-1,j)) * G%IdxT(i,j)
          elseif (valid_E) then ! Use a one-sided estimate from the east.
            sx = (S(i+1,j) - S(i,j)) * G%IdxCu(I,j)
          elseif (valid_W) then ! Use a one-sided estimate from the west.
            sx = (S(i,j) - S(i-1,j)) * G%IdxCu(I-1,j)
          endif
        endif

        ! Calculate the y-direction surface slope at tracer points.
        sy = 0.0
        valid_N = (ISS%hmask(i,j+1) == 1 .or. ISS%hmask(i,j+1) == 3)
        valid_S = (ISS%hmask(i,j-1) == 1 .or. ISS%hmask(i,j-1) == 3)
        if (CS%shelf_top_slope_bugs) then
          if (((j+j_off) == gjsc) .and. (.not. CS%reentrant_y)) then ! at south computational bdry
            if (valid_N) sy = (S(i,j+1)-S(i,j)) / G%dyT(i,j)
          elseif (((j+j_off) == gjec) .and. (.not. CS%reentrant_y)) then ! at north computational bdry
            if (valid_S) sy = (S(i,j)-S(i,j-1)) / G%dyT(i,j)
          elseif (valid_N .and. valid_S) then
            ! This is the usual interior point
            sy = (S(i,j+1) - S(i,j-1)) / (G%dyT(i,j) + G%dyT(i,j-1))
          elseif (valid_N) then
            sy = (S(i,j+1) - S(i,j)) / (G%dyT(i,j) + G%dyT(i,j+1))
          elseif (valid_S) then
            sy = (S(i,j) - S(i,j-1)) / (G%dyT(i,j) + G%dyT(i,j-1))
          endif
        else ! Correct the bugs in the version above.
          if (((j+j_off) == gjsc) .and. (.not. CS%reentrant_y)) then ! at south computational bdry
            if (valid_N) sy = (S(i,j+1) - S(i,j)) * G%IdyCv(i,J)
          elseif (((j+j_off) == gjec) .and. (.not. CS%reentrant_y)) then ! at north computational bdry
            if (valid_S) sy = (S(i,j) - S(i,j-1)) * G%IdyCv(i,J-1)
          elseif (valid_N .and. valid_S) then
            ! This is the usual interior point
            sy = 0.5*(S(i,j+1) - S(i,j-1)) * G%IdyT(i,j)
          elseif (valid_N) then ! Use a one-sided estimate from the north.
            sy = (S(i,j+1) - S(i,j)) * G%IdyCv(i,J)
          elseif (valid_S) then ! Use a one-sided estimate from the south.
            sy = (S(i,j) - S(i,j-1)) * G%IdyCv(i,J-1)
          endif
        endif

        if (CS%max_surface_slope>0) then
          scale = CS%max_surface_slope / max( sqrt((sx**2) + (sy**2)), CS%max_surface_slope )
          sx = scale*sx ; sy = scale*sy
        endif

        sx_e(i,j) = (-.25 * G%areaT(i,j)) * ((rho * grav) * (max(ISS%h_shelf(i,j),CS%min_h_shelf) * sx))
        sy_e(i,j) = (-.25 * G%areaT(i,j)) * ((rho * grav) * (max(ISS%h_shelf(i,j),CS%min_h_shelf) * sy))

        CS%sx_shelf(i,j) = sx ; CS%sy_shelf(i,j) = sy

        !Stress (Neumann) boundary conditions
        if (CS%ground_frac(i,j) == 1) then
          neumann_val = ((.5 * grav) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2 - &
                                        rhow * max(0.0, CS%bed_elev(i,j))**2))
        else
          neumann_val = (.5 * grav) * ((1-rho/rhow) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2))
        endif
        if ((CS%u_face_mask_bdry(I-1,j) == 2) .OR. &
          ((ISS%hmask(i-1,j) == 0 .OR. ISS%hmask(i-1,j) == 2) .AND. (CS%reentrant_x .OR. (i+i_off /= gisc)))) then
          ! left face of the cell is at a stress boundary
          ! the depth-integrated longitudinal stress is equal to the difference of depth-integrated
          ! pressure on either side of the face
          ! on the ice side, it is rho g h^2 / 2
          ! on the ocean side, it is rhow g (delta OD)^2 / 2
          ! OD can be zero under the ice; but it is ASSUMED on the ice-free side of the face, topography elevation
          !     is not above the base of the ice in the current cell

          ! Note the negative sign due to the direction of the normal vector
          taudx(I-1,J-1) = taudx(I-1,J-1) - .5 * G%dyCu(I-1,j) * neumann_val
          taudx(I-1,J) = taudx(I-1,J) - .5 * G%dyCu(I-1,j) * neumann_val
        endif

        if ((CS%u_face_mask_bdry(I,j) == 2) .OR. &
          ((ISS%hmask(i+1,j) == 0 .OR. ISS%hmask(i+1,j) == 2) .and. (CS%reentrant_x .OR. (i+i_off /= giec)))) then
          ! east face of the cell is at a stress boundary
          taudx(I,J-1) = taudx(I,J-1) + .5 * G%dyCu(I,j) * neumann_val
          taudx(I,J) = taudx(I,J) + .5 * G%dyCu(I,j) * neumann_val
        endif

        if ((CS%v_face_mask_bdry(i,J-1) == 2) .OR. &
          ((ISS%hmask(i,j-1) == 0 .OR. ISS%hmask(i,j-1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjsc)))) then
          ! south face of the cell is at a stress boundary
          taudy(I-1,J-1) = taudy(I-1,J-1) - .5 * G%dxCv(i,J-1) * neumann_val
          taudy(I,J-1) = taudy(I,J-1) - .5 * G%dxCv(i,J-1) * neumann_val
        endif

        if ((CS%v_face_mask_bdry(i,J) == 2) .OR. &
          ((ISS%hmask(i,j+1) == 0 .OR. ISS%hmask(i,j+1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjec)))) then
          ! north face of the cell is at a stress boundary
          taudy(I-1,J) = taudy(I-1,J) + .5 * G%dxCv(i,J) * neumann_val
          taudy(I,J) = taudy(I,J) + .5 * G%dxCv(i,J) * neumann_val
        endif
      else ! This is not an ice-filled cell, so zero out the slopes here
        CS%sx_shelf(i,j) = 0.0 ; CS%sy_shelf(i,j) = 0.0
        sx_e(i,j) = 0.0
        sy_e(i,j) = 0.0
      endif
    enddo
  enddo

  do J=jsc-1,jec ; do I=isc-1,iec
    taudx(I,J) = taudx(I,J) + ((sx_e(i,j)+sx_e(i+1,j+1)) + (sx_e(i+1,j)+sx_e(i,j+1)))
    taudy(I,J) = taudy(I,J) + ((sy_e(i,j)+sy_e(i+1,j+1)) + (sy_e(i+1,j)+sy_e(i,j+1)))
  enddo ; enddo
end subroutine calc_shelf_driving_stress

subroutine CG_action(CS, uret, vret, u_shlf, v_shlf, Phi, Phisub, umask, vmask, hmask, H_node, &
                     ice_visc, bathyT, u_curr, v_curr, G, US, is, ie, js, je, dens_ratio, &
                     use_newton_in, h_shelf)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: uret !< The retarding stresses working at u-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                         intent(inout) :: vret !< The retarding stresses working at v-points [R L3 Z T-2 ~> kg m s-2].
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                         intent(in)   :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(:,:,:,:,:,:), &
                         intent(in)    :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: u_shlf  !< The zonal ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: v_shlf  !< The meridional ice shelf velocity at vertices [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: umask !< A coded mask indicating the nature of the
                                             !! zonal flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: vmask !< A coded mask indicating the nature of the
                                             !! meridional flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G),CS%visc_qps), &
                         intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                               !! flow law [R L4 Z T-1 ~> kg m2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: bathyT !< The depth of ocean bathymetry at tracer points
                                                 !! relative to sea-level [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: u_curr  !< Frozen current iterate u^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(in)    :: v_curr  !< Frozen current iterate v^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]

  real,                  intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater, nondimensional
  type(unit_scale_type), intent(in)    :: US  !< A structure containing unit conversion factors
  integer,               intent(in)    :: is  !< The starting i-index to work on
  integer,               intent(in)    :: ie  !< The ending i-index to work on
  integer,               intent(in)    :: js  !< The starting j-index to work on
  integer,               intent(in)    :: je  !< The ending j-index to work on
  logical, optional,     intent(in)    :: use_newton_in !< If present, overrides CS%doing_newton for Newton correction
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         optional, intent(in) :: h_shelf !< Cell-averaged ice thickness for DG mode [Z ~> m]

! the linear action of the matrix on (u,v) with bilinear finite elements
! as of now everything is passed in so no grid pointers or anything of the sort have to be dereferenced,
! but this may change pursuant to conversations with others
!
! is & ie are the cells over which the iteration is done; this may change between calls to this subroutine
!     in order to make less frequent halo updates

! the linear action of the matrix on (u,v) with bilinear finite elements
! Phi has the form
! Phi(k,q,i,j) - applies to cell i,j

    !  3 - 4
    !  |   |
    !  1 - 2

! Phi(2*k-1,q,i,j) gives d(Phi_k)/dx at quadrature point q
! Phi(2*k,q,i,j) gives d(Phi_k)/dy at quadrature point q
! Phi_k is equal to 1 at vertex k, and 0 at vertex l /= k, and bilinear

  real :: ux, uy, vx, vy ! Components of velocity shears or divergence [T-1 ~> s-1]
  real :: uq, vq  ! Interpolated direction-vector δu at quadrature point [L T-1 ~> m s-1]
  real :: strx_n, stry_n, strsh_n, dstrain_n  ! Newton viscosity correction variables [T-1 ~> s-1], [T-2 ~> s-2]
  real :: u_curr_qp, v_curr_qp  ! Current iterate u^k at quadrature point [L T-1 ~> m s-1]
  real :: unorm2_qp  ! Regularized squared speed of u^k at quadrature point [L2 T-2 ~> m2 s-2]
  real :: basal_coef_qp  ! Picard basal friction coefficient at quadrature point [R L2 Z T-1 ~> kg s-1]
  real :: drag_newt_qp   ! Newton basal drag coefficient at quadrature point [R Z T-1 ~> kg m-2 s-1]
  real :: inner_dot_qp   ! u^k_qp · δu_qp inner product for Newton basal drag [L2 T-2 ~> m2 s-2]
  real :: coef_prefactor_e  ! Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real :: eps_vel2_e     ! Velocity regularization squared for current element [L2 T-2 ~> m2 s-2]
  real :: min_trac_e     ! min_basal_traction * areaT for current element [R L2 Z T-1 ~> kg s-1]
  real :: fB_e           ! Pre-computed Coulomb fB for element; 0 for Weertman [(T L-1)^CF_PostPeak]
  real :: jac_wt  ! Per-quadrature-point metric correction |J_q|/areaT [nondim]
  real :: h_gp           ! DG-evaluated ice thickness at Gauss point [Z ~> m]
  real :: bed_gp         ! Bilinear bed elevation at Gauss point [Z ~> m]
  real :: fB_local       ! Coulomb fB at quadrature point (DG mode) [(T L-1)^CF_PostPeak]
  real :: rho_oi_ratio   ! density_ocean / density_ice [nondim]
  real :: rho_ice_g_LtoZ ! US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  logical :: do_DG       ! Local flag for DG basal friction mode
  logical :: grounded_qp ! Whether this quadrature point is grounded (for DG per-qp check)
  real, dimension(:,:,:,:), pointer :: hgate ! Thickness field for the flotation
                         ! test: h_flot under DG_GL_GATE_CONTINUOUS, else h_nodal [Z ~> m]
  integer :: iq, jq, iphi, jphi, i, j, ilq, jlq, Itgt, Jtgt, qp, qpv
  logical :: visc_qp4
  logical :: use_newton  ! Whether to apply Newton tangent stiffness corrections
  logical :: do_newton_visc  ! Whether to apply viscosity-related Newton tangent stiffness corrections
  real, dimension(2) :: xquad  ! Nondimensional quadrature ratios [nondim]
  real, dimension(2,2) :: Usub, Vsub  ! Subgrid nodal contributions to basal traction [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2) :: Hcell   ! Ice shelf thickness at nodal (corner) points [Z ~> m]
  real, dimension(2,2,4) :: uret_qp, vret_qp                ! Temporary arrays in [R Z L3 T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: uret_b, vret_b  ! Temporary arrays in [R Z L3 T-2 ~> kg m s-2]

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  if (CS%visc_qps == 4) then
    visc_qp4=.true.
  else
    visc_qp4=.false.
    qpv = 1
  endif

  use_newton = CS%doing_newton
  if (present(use_newton_in)) use_newton = use_newton_in
  do_newton_visc = use_newton .and. trim(CS%ice_viscosity_compute) == "MODEL"

  do_DG = CS%use_DG_thickness .and. present(h_shelf)
  if (do_DG) then
    rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
    rho_ice_g_LtoZ = US%L_to_Z * CS%density_ice * CS%g_Earth
    hgate => CS%h_nodal
    if (CS%dg_gl_gate_continuous) hgate => CS%h_flot
  endif

  uret(:,:) = 0.0 ; vret(:,:) = 0.0
  uret_b(:,:,:) = 0.0 ; vret_b(:,:,:) = 0.0

  do j=js,je ; do i=is,ie ; if (hmask(i,j) == 1 .or. hmask(i,j)==3) then

    uret_qp(:,:,:) = 0.0 ; vret_qp(:,:,:) = 0.0

      ! Pre-computed element-level basal friction quantities (updated each outer Newton iteration
      ! by calc_shelf_basal_prefactors; avoids O(N_cg) recomputation of expensive prefactors).
      coef_prefactor_e = CS%coef_prefactor(i,j)
      eps_vel2_e = CS%eps_glen_min**2 * ((G%dxT(i,j)**2) + (G%dyT(i,j)**2))
      min_trac_e = CS%min_basal_traction * G%areaT(i,j)
      fB_e = CS%fB_elem(i,j)  ! 0 for Weertman; non-zero for Coulomb

      do iq=1,2 ; do jq=1,2

        qp = 2*(jq-1)+iq !current quad point

        uq = ((u_shlf(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
              (u_shlf(I,J) * (xquad(iq) * xquad(jq)))) + &
             ((u_shlf(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
              (u_shlf(I-1,J) * (xquad(3-iq) * xquad(jq))))

        vq = ((v_shlf(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
              (v_shlf(I,J) * (xquad(iq) * xquad(jq)))) + &
             ((v_shlf(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
              (v_shlf(I-1,J) * (xquad(3-iq) * xquad(jq))))

        ux = ((u_shlf(I-1,J-1) * Phi(1,qp,i,j)) + &
              (u_shlf(I,J) * Phi(7,qp,i,j))) + &
             ((u_shlf(I,J-1) * Phi(3,qp,i,j)) + &
              (u_shlf(I-1,J) * Phi(5,qp,i,j)))

        vx = ((v_shlf(I-1,J-1) * Phi(1,qp,i,j)) + &
              (v_shlf(I,J) * Phi(7,qp,i,j))) + &
             ((v_shlf(I,J-1) * Phi(3,qp,i,j)) + &
              (v_shlf(I-1,J) * Phi(5,qp,i,j)))

        uy = ((u_shlf(I-1,J-1) * Phi(2,qp,i,j)) + &
              (u_shlf(I,J) * Phi(8,qp,i,j))) + &
             ((u_shlf(I,J-1) * Phi(4,qp,i,j)) + &
              (u_shlf(I-1,J) * Phi(6,qp,i,j)))

        vy = ((v_shlf(I-1,J-1) * Phi(2,qp,i,j)) + &
              (v_shlf(I,J) * Phi(8,qp,i,j))) + &
             ((v_shlf(I,J-1) * Phi(4,qp,i,j)) + &
              (v_shlf(I-1,J) * Phi(6,qp,i,j)))

        if (visc_qp4) qpv = qp !current quad point for viscosity

        ! Newton correction: compute dstrain scalar once per quadrature point
        if (do_newton_visc) then
          strx_n = CS%newton_str_ux(i,j,qpv)
          stry_n = CS%newton_str_vy(i,j,qpv)
          strsh_n = CS%newton_str_sh(i,j,qpv)
          dstrain_n = (((2.*strx_n + stry_n)*ux) + ((2.*stry_n + strx_n)*vy)) + &
                      (strsh_n * (uy + vx) * 0.5)
        endif

        ! Basal friction and Newton Jacobian evaluated at this quadrature point (fully grounded cells only).
        ! Evaluating at quadrature points rather than cell-averaged ensures the Newton correction is the
        ! exact Jacobian of the Picard residual, enabling quadratic convergence for all friction exponents.
        grounded_qp = merge(CS%ground_frac(i,j) >= 1.0, CS%ground_frac(i,j) > 0.0, CS%GL_regularize)
        if (grounded_qp) then
          ! DG mode: per-Gauss-point grounding check and fB computation. h_gp is used
          ! only as a flotation measure (gate + effective pressure), so it is read from
          ! hgate (h_flot under DG_GL_GATE_CONTINUOUS).
          if (do_DG) then
            h_gp = ((hgate(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                    (hgate(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
                   ((hgate(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                    (hgate(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
            h_gp = max(h_gp, CS%min_h_shelf)
            bed_gp = ((CS%bed_node(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                      (CS%bed_node(I,J)     * (xquad(iq)   * xquad(jq))))  + &
                     ((CS%bed_node(I,J-1)   * (xquad(iq)   * xquad(3-jq))) + &
                      (CS%bed_node(I-1,J)   * (xquad(3-iq) * xquad(jq))))
            grounded_qp = (dens_ratio * h_gp - bed_gp > 0)
            if (grounded_qp .and. CS%CoulombFriction) then
              fB_local = compute_fB_local(h_gp, bed_gp, rho_oi_ratio, rho_ice_g_LtoZ, &
                  CS%C_basal_friction(i,j), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
                  CS%CF_PostPeak, CS%n_basal_fric)
            else
              fB_local = 0.0
            endif
          else
            fB_local = fB_e
          endif
        endif

        if (grounded_qp) then
          u_curr_qp = ((u_curr(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                       (u_curr(I,J) * (xquad(iq) * xquad(jq)))) + &
                      ((u_curr(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
                       (u_curr(I-1,J) * (xquad(3-iq) * xquad(jq))))
          v_curr_qp = ((v_curr(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                       (v_curr(I,J) * (xquad(iq) * xquad(jq)))) + &
                      ((v_curr(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
                       (v_curr(I-1,J) * (xquad(3-iq) * xquad(jq))))
          unorm2_qp = ((u_curr_qp**2) + (v_curr_qp**2)) + eps_vel2_e
          call compute_basal_coef(unorm2_qp, coef_prefactor_e, min_trac_e, fB_local, &
              CS%n_basal_fric, CS%CoulombFriction, CS%CF_PostPeak, US%L_T_to_m_s, use_newton, &
              basal_coef_qp, drag_newt_qp)
          ! Apply ground fraction scaling (replaces external scaling of basal_traction).
          ! Under GL_regularize, GL cells (0 < ground_frac < 1) get sub-grid-aware basal
          ! handling via the GL branch below, so the cell-level scaling must collapse to 1.0
          ! there to avoid double-counting (matches pre-refactor behavior when ground_frac
          ! was forced to 1.0 in GL cells).
          basal_coef_qp = basal_coef_qp * merge(1.0, CS%ground_frac(i,j), CS%GL_regularize)
          if (use_newton) then
            drag_newt_qp = drag_newt_qp * merge(1.0, CS%ground_frac(i,j), CS%GL_regularize)
            ! Inner product u^k_qp . delta_u_qp for the Newton correction.
            inner_dot_qp = (u_curr_qp * uq) + (v_curr_qp * vq)
          endif
        endif

        ! Ratio |J_q|/areaT corrects the uniform-area weight baked into ice_visc for
        ! non-rectangular elements where opposite cell edges have unequal lengths.
        jac_wt = CS%Jac(qp,i,j) * G%IareaT(i,j)

        do jphi=1,2 ; Jtgt = J-2+jphi ; do iphi=1,2 ; Itgt = I-2+iphi
          if (umask(Itgt,Jtgt) == 1) uret_qp(iphi,jphi,qp) = jac_wt * ice_visc(i,j,qpv) * &
            (((4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((uy+vx) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))
          if (vmask(Itgt,Jtgt) == 1) vret_qp(iphi,jphi,qp) = jac_wt * ice_visc(i,j,qpv) * &
            (((uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))

          ! Newton viscosity tangent stiffness: 2*(dη/dε_e^2) * (g·δε) * (g·φ_m).
          if (do_newton_visc) then
            if (umask(Itgt,Jtgt) == 1) uret_qp(iphi,jphi,qp) = uret_qp(iphi,jphi,qp) + &
              jac_wt * CS%newton_visc_factor(i,j,qpv) * dstrain_n * &
              (((2.*strx_n + stry_n) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
               (strsh_n * 0.5 * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))
            if (vmask(Itgt,Jtgt) == 1) vret_qp(iphi,jphi,qp) = vret_qp(iphi,jphi,qp) + &
              jac_wt * CS%newton_visc_factor(i,j,qpv) * dstrain_n * &
              ((strsh_n * 0.5 * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
               ((2.*stry_n + strx_n) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))
          endif

          if (grounded_qp) then
            ilq = 1 ; if (iq == iphi) ilq = 2
            jlq = 1 ; if (jq == jphi) jlq = 2
            ! Picard basal drag: C*|u^k|^(m-1) * δu evaluated at quadrature point, weighted by φ_m
            if (umask(Itgt,Jtgt) == 1) uret_qp(iphi,jphi,qp) = uret_qp(iphi,jphi,qp) + &
              (jac_wt * (basal_coef_qp * uq) * (xquad(ilq) * xquad(jlq)))
            if (vmask(Itgt,Jtgt) == 1) vret_qp(iphi,jphi,qp) = vret_qp(iphi,jphi,qp) + &
              (jac_wt * (basal_coef_qp * vq) * (xquad(ilq) * xquad(jlq)))
            ! Newton basal drag: pointwise Jacobian of the Picard residual.
            ! Tangent stiffness = basal_coef_qp*I + drag_newt_qp * u^k_qp ⊗ u^k_qp
            if (use_newton) then
              if (umask(Itgt,Jtgt) == 1) uret_qp(iphi,jphi,qp) = uret_qp(iphi,jphi,qp) + &
                jac_wt * drag_newt_qp * u_curr_qp * inner_dot_qp * (xquad(ilq) * xquad(jlq))
              if (vmask(Itgt,Jtgt) == 1) vret_qp(iphi,jphi,qp) = vret_qp(iphi,jphi,qp) + &
                jac_wt * drag_newt_qp * v_curr_qp * inner_dot_qp * (xquad(ilq) * xquad(jlq))
            endif
          endif
        enddo ; enddo
      enddo ; enddo

      !element contribution to SW node (node 1, which sees the current element as element 4)
      uret_b(I-1,J-1,4) = 0.25*((uret_qp(1,1,1)+uret_qp(1,1,4))+(uret_qp(1,1,2)+uret_qp(1,1,3)))
      vret_b(I-1,J-1,4) = 0.25*((vret_qp(1,1,1)+vret_qp(1,1,4))+(vret_qp(1,1,2)+vret_qp(1,1,3)))

      !element contribution to NW node (node 3, which sees the current element as element 2)
      uret_b(I-1,J  ,2) = 0.25*((uret_qp(1,2,1)+uret_qp(1,2,4))+(uret_qp(1,2,2)+uret_qp(1,2,3)))
      vret_b(I-1,J  ,2) = 0.25*((vret_qp(1,2,1)+vret_qp(1,2,4))+(vret_qp(1,2,2)+vret_qp(1,2,3)))

      !element contribution to SE node (node 2, which sees the current element as element 3)
      uret_b(I  ,J-1,3) = 0.25*((uret_qp(2,1,1)+uret_qp(2,1,4))+(uret_qp(2,1,2)+uret_qp(2,1,3)))
      vret_b(I  ,J-1,3) = 0.25*((vret_qp(2,1,1)+vret_qp(2,1,4))+(vret_qp(2,1,2)+vret_qp(2,1,3)))

      !element contribution to NE node (node 4, which sees the current element as element 1)
      uret_b(I  ,J  ,1) = 0.25*((uret_qp(2,2,1)+uret_qp(2,2,4))+(uret_qp(2,2,2)+uret_qp(2,2,3)))
      vret_b(I  ,J  ,1) = 0.25*((vret_qp(2,2,1)+vret_qp(2,2,4))+(vret_qp(2,2,2)+vret_qp(2,2,3)))

      if (CS%GL_regularize .and. CS%ground_frac(i,j) > 0.0 .and. CS%ground_frac(i,j) < 1.0) then
        ! Subgrid grounding-line: evaluate basal friction at each grounded sub-quadrature point.
        ! Picard and Newton Jacobian are both computed inside CG_action_subgrid_basal.
        Hcell(:,:) = H_node(I-1:I,J-1:J)
        if (do_DG) then
          ! h_nodal_cell is used inside only as a flotation measure (sub-qp gate +
          ! effective pressure), so the gate field is passed (h_flot under
          ! DG_GL_GATE_CONTINUOUS).
          call CG_action_subgrid_basal(CS, G, US, Phisub, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
              bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
              use_DG=.true., h_shelf_cell=h_shelf(i,j), &
              h_nodal_cell=hgate(i,j,:,:), &
              bed_corners=CS%bed_node(I-1:I,J-1:J))
        else
          call CG_action_subgrid_basal(CS, G, US, Phisub, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
              bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j))
        endif
        if (umask(I-1,J-1) == 1) uret_b(I-1,J-1,4) = uret_b(I-1,J-1,4) + Usub(1,1)
        if (umask(I-1,J  ) == 1) uret_b(I-1,J  ,2) = uret_b(I-1,J  ,2) + Usub(1,2)
        if (umask(I  ,J-1) == 1) uret_b(I  ,J-1,3) = uret_b(I  ,J-1,3) + Usub(2,1)
        if (umask(I  ,J  ) == 1) uret_b(I  ,J  ,1) = uret_b(I  ,J  ,1) + Usub(2,2)
        if (vmask(I-1,J-1) == 1) vret_b(I-1,J-1,4) = vret_b(I-1,J-1,4) + Vsub(1,1)
        if (vmask(I-1,J  ) == 1) vret_b(I-1,J  ,2) = vret_b(I-1,J  ,2) + Vsub(1,2)
        if (vmask(I  ,J-1) == 1) vret_b(I  ,J-1,3) = vret_b(I  ,J-1,3) + Vsub(2,1)
        if (vmask(I  ,J  ) == 1) vret_b(I  ,J  ,1) = vret_b(I  ,J  ,1) + Vsub(2,2)
      endif
  endif ; enddo ; enddo

  do J=js-1,je ; do I=is-1,ie
    uret(I,J) = (uret_b(I,J,1)+uret_b(I,J,4)) + (uret_b(I,J,2)+uret_b(I,J,3))
    vret(I,J) = (vret_b(I,J,1)+vret_b(I,J,4)) + (vret_b(I,J,2)+vret_b(I,J,3))
  enddo ; enddo

end subroutine CG_action

!> Compute subgrid grounding-line basal traction nodal contributions for a CG action.
!! Evaluates basal friction (Picard and Newton Jacobian) at each grounded sub-quadrature point.
!! The sub-qp flotation test accounts for partial grounding; no external ground_frac scaling needed.
subroutine CG_action_subgrid_basal(CS, G, US, Phisub, H, U_curr, V_curr, U_delta, V_delta, &
                                   bathyT, dens_ratio, i_elem, j_elem, fB_e, use_newton, Ucontr, Vcontr, &
                                   dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT, &
                                   use_DG, h_shelf_cell, h_nodal_cell, bed_corners)
  type(ice_shelf_dyn_CS), intent(in) :: CS      !< Ice shelf control structure
  type(ocean_grid_type),  intent(in) :: G       !< The grid structure
  type(unit_scale_type),  intent(in) :: US      !< Unit conversion factors
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Sub-grid quadrature weights [nondim]
  real, dimension(2,2),   intent(in) :: H       !< Ice thickness at element corners [Z ~> m]
  real, dimension(2,2),   intent(in) :: U_curr  !< Frozen u^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_curr  !< Frozen v^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: U_delta !< Search direction δu at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_delta !< Search direction δv at element corners [L T-1 ~> m s-1]
  real,                   intent(in) :: bathyT  !< Ocean bathymetry depth at tracer point [Z ~> m]
  real,                   intent(in) :: dens_ratio !< Ice density / water density [nondim]
  integer,                intent(in) :: i_elem  !< Tracer-grid i-index of the element
  integer,                intent(in) :: j_elem  !< Tracer-grid j-index of the element
  real,                   intent(in) :: fB_e    !< Element Coulomb parameter fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  logical,                intent(in) :: use_newton !< If true, include Newton basal drag correction
  real, dimension(2,2),   intent(out) :: Ucontr !< Nodal u-contributions with friction applied [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2),   intent(out) :: Vcontr !< Nodal v-contributions with friction applied [R L3 Z T-2 ~> kg m s-2]
  real,                   intent(in) :: dxCv_S !< The cell width at the southern (v-point) edge [L ~> m]
  real,                   intent(in) :: dxCv_N !< The cell width at the northern (v-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_W !< The cell height at the western (u-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_E !< The cell height at the eastern (u-point) edge [L ~> m]
  real,                   intent(in) :: IareaT !< The inverse of the cell area at the tracer point [L-2 ~> m-2]
  logical,       optional, intent(in) :: use_DG       !< If true, use DG thickness and bed_node [nondim]
  real,          optional, intent(in) :: h_shelf_cell  !< Cell-averaged ice thickness [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: h_nodal_cell !< Q1 thickness at the 4 cell corners, used
                                          !! only as a flotation measure (sub-qp gate + effective
                                          !! pressure); the caller passes h_flot here under
                                          !! DG_GL_GATE_CONTINUOUS [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: bed_corners !< Bed elevation at element corners [Z ~> m]

  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: Ucontr_sub, Vcontr_sub
  real, dimension(2,2,2,2) :: U_qp_nd, V_qp_nd  ! Per-qp nodal contributions (qx,qy,m,n)
                                                ! accumulated then pair-summed for rotation invariance
  real :: hloc          ! Local sub-cell ice thickness [Z ~> m]
  real :: bed_sub       ! Bed elevation at sub-quadrature point [Z ~> m]
  real :: u_curr_loc    ! Frozen u^k interpolated to sub-qp [L T-1 ~> m s-1]
  real :: v_curr_loc    ! Frozen v^k interpolated to sub-qp [L T-1 ~> m s-1]
  real :: u_delta_loc   ! Search direction δu interpolated to sub-qp [L T-1 ~> m s-1]
  real :: v_delta_loc   ! Search direction δv interpolated to sub-qp [L T-1 ~> m s-1]
  real :: unorm2_loc    ! Regularized |u^k|^2 at sub-qp [L2 T-2 ~> m2 s-2]
  real :: basal_coef_loc ! Picard friction coefficient at sub-qp [R L2 Z T-1 ~> kg s-1]
  real :: drag_newt_loc  ! Newton drag coefficient at sub-qp [R Z T ~> kg m-2 s]
  real :: inner_dot_loc  ! u^k · δu inner product at sub-qp [L2 T-2 ~> m2 s-2]
  real :: phi_mn         ! Basis function value at sub-qp [nondim]
  real :: contrib        ! Quadrature weight contribution [nondim]
  real :: coef_prefactor ! Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real :: min_trac_area  ! Minimum area-integrated traction floor [R L2 Z T-1 ~> kg s-1]
  real :: eps_vel2       ! Velocity regularization squared [L2 T-2 ~> m2 s-2]
  real :: jac_sub_wt ! Per-sub-cell-QP metric correction |J_sub|/areaT [nondim]
  real :: a, d      ! Interpolated cell-edge spacings at the sub-cell QP [L ~> m]
  real :: subarea        ! Fractional sub-cell area [nondim]
  real :: fB_local       ! Coulomb fB at sub-qp (DG mode) [(T L-1)^CF_PostPeak]
  real :: rho_oi_ratio   ! density_ocean / density_ice [nondim]
  real :: rho_ice_g_LtoZ ! US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  real :: xi_sub, eta_sub ! DG reference coords at sub-qp ([-0.5,0.5]) [nondim]
  logical :: do_DG       ! Local flag for DG mode
  integer :: nsub, i, j, qx, qy, m, n

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  coef_prefactor = CS%coef_prefactor(i_elem,j_elem)
  min_trac_area  = CS%min_basal_traction * G%areaT(i_elem,j_elem)
  eps_vel2 = CS%eps_glen_min**2 * ((G%dxT(i_elem,j_elem)**2) + (G%dyT(i_elem,j_elem)**2))

  do_DG = .false.
  if (present(use_DG)) do_DG = use_DG
  if (do_DG) then
    rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
    rho_ice_g_LtoZ = US%L_to_Z * CS%density_ice * CS%g_Earth
  endif

  Ucontr_sub(:,:,:,:) = 0.0 ; Vcontr_sub(:,:,:,:) = 0.0

  do j=1,nsub ; do i=1,nsub
    U_qp_nd(:,:,:,:) = 0.0 ; V_qp_nd(:,:,:,:) = 0.0
    do qy=1,2 ; do qx=1,2
      if (do_DG) then
        ! xi_sub = a_right(qx,i) - 0.5; marginal sum of Phisub over the l index gives a_right(qx,i).
        xi_sub  = (Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)) - 0.5
        eta_sub = (Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)) - 0.5
        ! Nodal Q1 evaluation at the sub-QP via the Phisub corner-basis weights.
        hloc = ((Phisub(qx,qy,i,j,1,1)*h_nodal_cell(1,1)) + (Phisub(qx,qy,i,j,2,2)*h_nodal_cell(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*h_nodal_cell(1,2)) + (Phisub(qx,qy,i,j,2,1)*h_nodal_cell(2,1)))
        hloc = max(hloc, CS%min_h_shelf)
        bed_sub = ((Phisub(qx,qy,i,j,1,1)*bed_corners(1,1)) + (Phisub(qx,qy,i,j,2,2)*bed_corners(2,2))) + &
                  ((Phisub(qx,qy,i,j,1,2)*bed_corners(1,2)) + (Phisub(qx,qy,i,j,2,1)*bed_corners(2,1)))
      else
        ! Standard mode: bilinear H interpolation, cell-averaged bed
        hloc = ((Phisub(qx,qy,i,j,1,1)*H(1,1)) + (Phisub(qx,qy,i,j,2,2)*H(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*H(1,2)) + (Phisub(qx,qy,i,j,2,1)*H(2,1)))
        bed_sub = bathyT
      endif

      if (dens_ratio * hloc - bed_sub > 0) then  ! grounded sub-qp
        u_curr_loc  = (((Phisub(qx,qy,i,j,1,1)*U_curr(1,1))  + (Phisub(qx,qy,i,j,2,2)*U_curr(2,2)))  + &
                       ((Phisub(qx,qy,i,j,1,2)*U_curr(1,2))  + (Phisub(qx,qy,i,j,2,1)*U_curr(2,1))))
        v_curr_loc  = (((Phisub(qx,qy,i,j,1,1)*V_curr(1,1))  + (Phisub(qx,qy,i,j,2,2)*V_curr(2,2)))  + &
                       ((Phisub(qx,qy,i,j,1,2)*V_curr(1,2))  + (Phisub(qx,qy,i,j,2,1)*V_curr(2,1))))
        u_delta_loc = (((Phisub(qx,qy,i,j,1,1)*U_delta(1,1)) + (Phisub(qx,qy,i,j,2,2)*U_delta(2,2))) + &
                       ((Phisub(qx,qy,i,j,1,2)*U_delta(1,2)) + (Phisub(qx,qy,i,j,2,1)*U_delta(2,1))))
        v_delta_loc = (((Phisub(qx,qy,i,j,1,1)*V_delta(1,1)) + (Phisub(qx,qy,i,j,2,2)*V_delta(2,2))) + &
                       ((Phisub(qx,qy,i,j,1,2)*V_delta(1,2)) + (Phisub(qx,qy,i,j,2,1)*V_delta(2,1))))

        unorm2_loc = ((u_curr_loc**2) + (v_curr_loc**2)) + eps_vel2

        ! Compute Coulomb fB at this sub-qp when in DG mode
        if (do_DG .and. CS%CoulombFriction) then
          fB_local = compute_fB_local(hloc, bed_sub, rho_oi_ratio, rho_ice_g_LtoZ, &
              CS%C_basal_friction(i_elem,j_elem), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
              CS%CF_PostPeak, CS%n_basal_fric)
        else
          fB_local = fB_e
        endif

        call compute_basal_coef(unorm2_loc, coef_prefactor, min_trac_area, fB_local, &
            CS%n_basal_fric, CS%CoulombFriction, CS%CF_PostPeak, US%L_T_to_m_s, use_newton, &
            basal_coef_loc, drag_newt_loc)
        inner_dot_loc = (u_curr_loc * u_delta_loc) + (v_curr_loc * v_delta_loc)

        ! Interpolate cell-edge metrics to the sub-cell QP using the bilinear shape function values
        ! from bilinear_shape_functions_subgrid.  Marginal sums of Phisub give the interpolation
        ! weights: sum over k=1 nodes gives (1-y); k=2 gives y; l=1 gives (1-x); l=2 gives x.
        ! This is analogous to jac_wt = CS%Jac(qp,i,j) * G%IareaT(i,j) in the regular routines.
        a = (dxCv_S * (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1))) + &  ! (1-y) * dxCv_S
            (dxCv_N * (Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)))      !  + y  * dxCv_N
        d = (dyCu_W * (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2))) + &  ! (1-x) * dyCu_W
            (dyCu_E * (Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)))      !  + x  * dyCu_E
        jac_sub_wt = 0.25 * subarea * (a * d) * IareaT

        do n=1,2 ; do m=1,2
          phi_mn  = Phisub(qx,qy,i,j,m,n)
          contrib = jac_sub_wt * phi_mn
          ! Picard: friction matrix applied to search direction δu
          U_qp_nd(qx,qy,m,n) = contrib * (basal_coef_loc * u_delta_loc)
          V_qp_nd(qx,qy,m,n) = contrib * (basal_coef_loc * v_delta_loc)
          ! Newton: Jacobian d(tau_b_i)/d(u_j) = basal_coef*I + drag_newt*u^k_i*u^k_j
          if (use_newton) then
            U_qp_nd(qx,qy,m,n) = U_qp_nd(qx,qy,m,n) + (contrib * (drag_newt_loc * u_curr_loc * inner_dot_loc))
            V_qp_nd(qx,qy,m,n) = V_qp_nd(qx,qy,m,n) + (contrib * (drag_newt_loc * v_curr_loc * inner_dot_loc))
          endif
        enddo ; enddo
      endif
    enddo ; enddo

    do n=1,2 ; do m=1,2
      Ucontr_sub(i,j,m,n) = (U_qp_nd(1,1,m,n) + U_qp_nd(2,2,m,n)) + &
                            (U_qp_nd(1,2,m,n) + U_qp_nd(2,1,m,n))
      Vcontr_sub(i,j,m,n) = (V_qp_nd(1,1,m,n) + V_qp_nd(2,2,m,n)) + &
                            (V_qp_nd(1,2,m,n) + V_qp_nd(2,1,m,n))
    enddo ; enddo
  enddo ; enddo

  do n=1,2 ; do m=1,2
    call sum_square_matrix(Ucontr(m,n), Ucontr_sub(:,:,m,n), nsub)
    call sum_square_matrix(Vcontr(m,n), Vcontr_sub(:,:,m,n), nsub)
  enddo ; enddo

end subroutine CG_action_subgrid_basal

!> Compute the Picard basal friction coefficient and Newton drag coefficient at a
!! single quadrature point. Encapsulates the 3-path dispatch (linear Weertman / nonlinear
!! Weertman / Coulomb) so that CG_action, matrix_diagonal, and their subgrid equivalents
!! remain readable. The ground_frac scaling is NOT applied here; callers do it after the call.
subroutine compute_basal_coef(unorm2_qp, coef_prefactor, min_trac_area, fB_e, &
    n_basal_fric, CoulombFriction, CF_PostPeak, L_T_to_m_s, use_newton, &
    basal_coef, drag_newt)
  real,    intent(in)  :: unorm2_qp      !< Regularized |u^k|^2 > 0 at quadrature point [L2 T-2 ~> m2 s-2]
  real,    intent(in)  :: coef_prefactor !< Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real,    intent(in)  :: min_trac_area  !< Pre-computed min_basal_traction * areaT floor [R L2 Z T-1 ~> kg s-1]
  real,    intent(in)  :: fB_e           !< Element-level Coulomb fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  real,    intent(in)  :: n_basal_fric   !< Friction sliding exponent m [nondim]
  logical, intent(in)  :: CoulombFriction !< True if using Coulomb friction
  real,    intent(in)  :: CF_PostPeak    !< Coulomb post-peak exponent q [nondim]
  real,    intent(in)  :: L_T_to_m_s    !< Unit conversion factor from internal [L T-1] to [m s-1]
  logical, intent(in)  :: use_newton     !< If true, evaluate drag_newt; otherwise set to 0
  real,    intent(out) :: basal_coef     !< Picard friction coefficient at quadrature point [R L2 Z T-1 ~> kg s-1]
  real,    intent(out) :: drag_newt      !< Newton drag coefficient [R Z T ~> kg m-2 s]; 0 without Newton

  real :: unorm    ! |u^k| at quadrature point in physical units [m s-1]
  real :: raw_coef ! Pre-floor friction coefficient [R L2 Z T-1 ~> kg s-1]
  real :: fBuq     ! fB_e * |u^k|^q [nondim]

  if (n_basal_fric == 1.0 .and. .not. CoulombFriction) then
    ! Linear Weertman: coef is independent of |u|; sqrt and Newton correction not needed
    basal_coef = max(coef_prefactor, min_trac_area)
    drag_newt  = 0.0
  elseif (CoulombFriction) then
    ! Schoof/Gagliardini Coulomb friction
    unorm    = L_T_to_m_s * sqrt(unorm2_qp)
    fBuq     = fB_e * unorm**CF_PostPeak
    raw_coef = coef_prefactor * (unorm**(n_basal_fric-1.0)) / (1.0 + fBuq)**n_basal_fric
    if (raw_coef < min_trac_area) then
      basal_coef = min_trac_area  ;  drag_newt = 0.0
    else
      basal_coef = raw_coef
      if (use_newton) then
        drag_newt = (1.0/unorm2_qp) * raw_coef * &
            ((n_basal_fric-1.0) - n_basal_fric * CF_PostPeak * fBuq / (1.0 + fBuq))
      else
        drag_newt = 0.0
      endif
    endif
  else
    ! Nonlinear Weertman (m > 1)
    unorm    = L_T_to_m_s * sqrt(unorm2_qp)
    raw_coef = coef_prefactor * (unorm**(n_basal_fric-1.0))
    if (raw_coef < min_trac_area) then
      basal_coef = min_trac_area  ;  drag_newt = 0.0
    else
      basal_coef = raw_coef
      if (use_newton) then
        drag_newt = (n_basal_fric-1.0) / unorm2_qp * raw_coef
      else
        drag_newt = 0.0
      endif
    endif
  endif

end subroutine compute_basal_coef

!> Compute the Coulomb fB parameter at a single quadrature point from local (subgrid) ice
!! thickness and bed elevation. This replaces the cell-averaged fB_elem when use_DG_thickness
!! is active, allowing the effective pressure to vary within the cell.
pure real function compute_fB_local(h_local, bed_local, rho_oi_ratio, rho_ice_g_LtoZ, &
    C_basal, alpha_coulomb, CF_Max, CF_MinN, CF_PostPeak, n_basal_fric)
  real, intent(in) :: h_local       !< Ice thickness at quadrature point [Z ~> m]
  real, intent(in) :: bed_local     !< Bed elevation at quadrature point [Z ~> m]
  real, intent(in) :: rho_oi_ratio  !< density_ocean_avg / density_ice [nondim]
  real, intent(in) :: rho_ice_g_LtoZ !< US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  real, intent(in) :: C_basal       !< Basal friction coefficient for this cell [R L Z T-2 (s m-1)^n]
  real, intent(in) :: alpha_coulomb !< Coulomb prefactor [nondim]
  real, intent(in) :: CF_Max        !< Coulomb friction maximum coefficient [nondim]
  real, intent(in) :: CF_MinN       !< Minimum Coulomb effective pressure [R Z L T-2 ~> Pa]
  real, intent(in) :: CF_PostPeak   !< Coulomb post-peak exponent q [nondim]
  real, intent(in) :: n_basal_fric  !< Friction sliding exponent m [nondim]

  real :: Hf  ! Flotation thickness [Z ~> m]
  real :: fN  ! Effective pressure [R Z L T-2 ~> Pa]

  Hf = max(rho_oi_ratio * bed_local, 0.0)
  fN = max(rho_ice_g_LtoZ * (h_local - Hf), CF_MinN)
  compute_fB_local = alpha_coulomb * (C_basal / (CF_Max * fN))**(CF_PostPeak / n_basal_fric)
end function compute_fB_local

!! Returns the sum of the elements in a square matrix. This sum is bitwise identical even if the matrices are rotated.
subroutine sum_square_matrix(sum_out, mat_in, n)
  integer, intent(in) :: n !< The length and width of each matrix in mat_in
  real, dimension(n,n), intent(in) :: mat_in !< The n x n matrix whose elements will be summed
  real, intent(out) :: sum_out !< The sum of the elements of matrix mat_in
  integer :: s0, e0, s1, e1

  sum_out = 0.0

  s0 = 1 ; e0 = n

  !start by summing elements on outer edges of matrix
  do while (s0<e0)

    !corners
    sum_out = sum_out + ( (mat_in(s0,s0) + mat_in(e0,e0)) + (mat_in(e0,s0) + mat_in(s0,e0)) )

    s1 = s0+1 ; e1 = e0-1

    do while (s1<e1) !non-corners

      sum_out = sum_out + &
                ( ( (mat_in(s0,s1) + mat_in(s1,s0)) + (mat_in(e0,e1) + mat_in(e1,e0)) ) + &
                  ( (mat_in(e1,s0) + mat_in(e0,s1)) + (mat_in(s1,e0) + mat_in(s0,e1)) ) )

      s1 = s1+1 ; e1 = e1-1
    enddo

    !center element of an edge
    if (s1==e1) sum_out = sum_out + ( (mat_in(s1,s0) + mat_in(e1,e0)) + (mat_in(e0,e1) + mat_in(s0,s1)) )

    s0 = s0+1 ; e0 = e0-1 !next loop iteration using new edges that are one element inward of the current edges
  enddo

  !center element of entire matrix
  if (s0==e0) sum_out = sum_out + mat_in(s0,e0)

end subroutine sum_square_matrix

!> returns the diagonal entries of the matrix for a Jacobi preconditioning
subroutine matrix_diagonal(CS, G, US, H_node, ice_visc, u_curr, v_curr, &
                           hmask, dens_ratio, Phi, Phisub, u_diagonal, v_diagonal, h_shelf)

  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< The ice shelf thickness at nodal
                                                 !! (corner) points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G),CS%visc_qps), &
                          intent(in)    :: ice_visc !< A field related to the ice viscosity from Glen's
                                                !! flow law [R L4 Z T-1 ~> kg m2 s-1].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: u_curr  !< Frozen current iterate u^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: v_curr  !< Frozen current iterate v^k, used to evaluate basal friction
                                               !! at quadrature points [L T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                   intent(in)    :: dens_ratio !< The density of ice divided by the density
                                                     !! of seawater [nondim]
  real, dimension(8,4,SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1]
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Quadrature structure weights at subgridscale
                                            !! locations for finite element calculations [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: u_diagonal !< The diagonal elements of the u-velocity
                                            !! matrix from the left-hand side of the solver [R L2 Z T-1 ~> kg s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: v_diagonal  !< The diagonal elements of the v-velocity
                                            !! matrix from the left-hand side of the solver [R L2 Z T-1 ~> kg s-1]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         optional, intent(in) :: h_shelf !< Cell-averaged ice thickness for DG mode [Z ~> m]


! returns the diagonal entries of the matrix for a Jacobi preconditioning

  real :: ux, uy, vx, vy ! Interpolated weight gradients [L-1 ~> m-1]
  real :: jac_wt  ! Per-quadrature-point metric correction |J_q|/areaT [nondim]
  real :: strx_n, stry_n, strsh_n  ! Newton viscosity strain rates [T-1 ~> s-1]
  real :: dstrain_diag_u, dstrain_diag_v  ! Newton viscosity diagonal correction factors [T-1 L-1 ~> s-1 m-1]
  real :: phi_m_sq  ! Squared basis function value at quadrature point [nondim]
  real :: u_curr_qp, v_curr_qp  ! Current iterate u^k at quadrature point [L T-1 ~> m s-1]
  real :: unorm2_qp  ! Regularized squared speed of u^k at quadrature point [L2 T-2 ~> m2 s-2]
  real :: basal_coef_qp  ! Picard basal friction coefficient at quadrature point [R L2 Z T-1 ~> kg s-1]
  real :: drag_newt_qp   ! Newton basal drag coefficient at quadrature point [R Z T-1 ~> kg m-2 s-1]
  real :: coef_prefactor_e  ! Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real :: eps_vel2_e     ! Velocity regularization squared for current element [L2 T-2 ~> m2 s-2]
  real :: min_trac_e     ! min_basal_traction * areaT for current element [R L2 Z T-1 ~> kg s-1]
  real :: fB_e           ! Pre-computed Coulomb fB for element; 0 for Weertman [(T L-1)^CF_PostPeak]
  real :: h_gp           ! DG-evaluated ice thickness at Gauss point [Z ~> m]
  real :: bed_gp         ! Bilinear bed elevation at Gauss point [Z ~> m]
  real :: fB_local       ! Coulomb fB at quadrature point (DG mode) [(T L-1)^CF_PostPeak]
  real :: rho_oi_ratio   ! density_ocean / density_ice [nondim]
  real :: rho_ice_g_LtoZ ! US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  logical :: do_DG       ! Local flag for DG basal friction mode
  logical :: grounded_qp ! Whether this quadrature point is grounded
  real, dimension(:,:,:,:), pointer :: hgate ! Thickness field for the flotation
                         ! test: h_flot under DG_GL_GATE_CONTINUOUS, else h_nodal [Z ~> m]
  real, dimension(2)   :: xquad
  real, dimension(2,2) :: Hcell, u_diag_sub, v_diag_sub  ! Subgrid diagonal contributions [R L2 Z T-1 ~> kg s-1]
  real, dimension(2,2,4) :: u_diag_qp, v_diag_qp
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: u_diag_b, v_diag_b
  logical :: do_newton_visc  ! Whether to apply viscosity-related Newton tangent stiffness corrections
  logical :: visc_qp4
  integer :: i, j, isc, jsc, iec, jec, iphi, jphi, iq, jq, ilq, jlq, Itgt, Jtgt, qp, qpv

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))

  if (CS%visc_qps == 4) then
    visc_qp4=.true.
  else
    visc_qp4=.false.
    qpv = 1
  endif

  do_newton_visc = CS%doing_newton .and. trim(CS%ice_viscosity_compute) == "MODEL"

  do_DG = CS%use_DG_thickness .and. present(h_shelf)
  if (do_DG) then
    rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
    rho_ice_g_LtoZ = US%L_to_Z * CS%density_ice * CS%g_Earth
    hgate => CS%h_nodal
    if (CS%dg_gl_gate_continuous) hgate => CS%h_flot
  endif

  u_diag_b(:,:,:)=0.0
  v_diag_b(:,:,:)=0.0

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1 ; if (hmask(i,j) == 1 .or. hmask(i,j)==3) then

    ! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
    ! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j

    u_diag_qp(:,:,:) = 0.0 ; v_diag_qp(:,:,:) = 0.0

      ! Pre-computed element-level basal friction quantities (updated each outer iteration).
      coef_prefactor_e = CS%coef_prefactor(i,j)
      eps_vel2_e = CS%eps_glen_min**2 * ((G%dxT(i,j)**2) + (G%dyT(i,j)**2))
      min_trac_e = CS%min_basal_traction * G%areaT(i,j)
      fB_e = CS%fB_elem(i,j)  ! 0 for Weertman; non-zero for Coulomb

    do iq=1,2 ; do jq=1,2

      qp = 2*(jq-1)+iq !current quad point
      if (visc_qp4) qpv = qp !current quad point for viscosity

      ! Ratio |J_q|/areaT corrects the uniform-area weight baked into ice_visc for
      ! non-rectangular elements where opposite cell edges have unequal lengths.
      jac_wt = CS%Jac(qp,i,j) * G%IareaT(i,j)

      ! Pre-compute Newton strain data for this QP (for viscosity diagonal correction)
      if (do_newton_visc) then
        strx_n = CS%newton_str_ux(i,j,qpv)
        stry_n = CS%newton_str_vy(i,j,qpv)
        strsh_n = CS%newton_str_sh(i,j,qpv)
      endif

      ! Basal friction coefficients at this quadrature point (fully grounded cells only)
      grounded_qp = merge(CS%ground_frac(i,j) >= 1.0, CS%ground_frac(i,j) > 0.0, CS%GL_regularize)
      if (grounded_qp) then
        if (do_DG) then
          ! h_gp is used only as a flotation measure (gate + effective pressure), so it
          ! is read from hgate (h_flot under DG_GL_GATE_CONTINUOUS).
          h_gp = ((hgate(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                  (hgate(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
                 ((hgate(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                  (hgate(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
          h_gp = max(h_gp, CS%min_h_shelf)
          ! Bed at the QP with the same rotation-paired bilinear pattern as
          ! u_curr_qp below, for symmetric rotation-cancel structure.
          bed_gp = ((CS%bed_node(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                    (CS%bed_node(I,J)     * (xquad(iq)   * xquad(jq))))  + &
                   ((CS%bed_node(I,J-1)   * (xquad(iq)   * xquad(3-jq))) + &
                    (CS%bed_node(I-1,J)   * (xquad(3-iq) * xquad(jq))))
          grounded_qp = (dens_ratio * h_gp - bed_gp > 0)
          if (grounded_qp .and. CS%CoulombFriction) then
            fB_local = compute_fB_local(h_gp, bed_gp, rho_oi_ratio, rho_ice_g_LtoZ, &
                CS%C_basal_friction(i,j), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
                CS%CF_PostPeak, CS%n_basal_fric)
          else
            fB_local = 0.0
          endif
        else
          fB_local = fB_e
        endif
      endif

      if (grounded_qp) then
        u_curr_qp = ((u_curr(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                     (u_curr(I,J) * (xquad(iq) * xquad(jq)))) + &
                    ((u_curr(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
                     (u_curr(I-1,J) * (xquad(3-iq) * xquad(jq))))
        v_curr_qp = ((v_curr(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                     (v_curr(I,J) * (xquad(iq) * xquad(jq)))) + &
                    ((v_curr(I,J-1) * (xquad(iq) * xquad(3-jq))) + &
                     (v_curr(I-1,J) * (xquad(3-iq) * xquad(jq))))
        unorm2_qp = ((u_curr_qp**2) + (v_curr_qp**2)) + eps_vel2_e
        call compute_basal_coef(unorm2_qp, coef_prefactor_e, min_trac_e, fB_local, &
            CS%n_basal_fric, CS%CoulombFriction, CS%CF_PostPeak, US%L_T_to_m_s, .true., &
            basal_coef_qp, drag_newt_qp)
        basal_coef_qp = basal_coef_qp * merge(1.0, CS%ground_frac(i,j), CS%GL_regularize)
        drag_newt_qp  = drag_newt_qp  * merge(1.0, CS%ground_frac(i,j), CS%GL_regularize)
      endif

      do jphi=1,2 ; Jtgt = J-2+jphi ; do iphi=1,2 ; Itgt = I-2+iphi

        ilq = 1 ; if (iq == iphi) ilq = 2
        jlq = 1 ; if (jq == jphi) jlq = 2
        phi_m_sq = (xquad(ilq) * xquad(jlq))**2

        if (CS%umask(Itgt,Jtgt) == 1) then

          ux = Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)
          uy = Phi(2*(2*(jphi-1)+iphi),qp,i,j)
          vx = 0.
          vy = 0.

          u_diag_qp(iphi,jphi,qp) = jac_wt * &
            ice_visc(i,j,qpv) * (((4*ux+2*vy) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((uy+vx) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))

          ! Newton viscosity diagonal correction: newton_visc_factor * (g . grad_phi_m_u)^2
          ! where grad_phi_m_u = [(2*strx+stry)*Phi_xm + strsh/2*Phi_ym] for u-DOF at node m
          if (do_newton_visc) then
            dstrain_diag_u = ((2.*strx_n + stry_n) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
                             (strsh_n * 0.5 * Phi(2*(2*(jphi-1)+iphi),qp,i,j))
            u_diag_qp(iphi,jphi,qp) = u_diag_qp(iphi,jphi,qp) + &
              jac_wt * CS%newton_visc_factor(i,j,qpv) * dstrain_diag_u**2
          endif

          if (grounded_qp) then
            u_diag_qp(iphi,jphi,qp) = u_diag_qp(iphi,jphi,qp) + jac_wt * basal_coef_qp * phi_m_sq
            if (CS%doing_newton) &
              u_diag_qp(iphi,jphi,qp) = u_diag_qp(iphi,jphi,qp) + jac_wt * drag_newt_qp * u_curr_qp**2 * phi_m_sq
          endif
        endif

        if (CS%vmask(Itgt,Jtgt) == 1) then

          vx = Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)
          vy = Phi(2*(2*(jphi-1)+iphi),qp,i,j)
          ux = 0.
          uy = 0.

          v_diag_qp(iphi,jphi,qp) = jac_wt *  &
            ice_visc(i,j,qpv) * (((uy+vx) * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
            ((4*vy+2*ux) * Phi(2*(2*(jphi-1)+iphi),qp,i,j)))

          ! Newton viscosity diagonal correction for v-DOF: uses [strsh/2*Phi_xm + (2*stry+strx)*Phi_ym]
          if (do_newton_visc) then
            dstrain_diag_v = (strsh_n * 0.5 * Phi(2*(2*(jphi-1)+iphi)-1,qp,i,j)) + &
                             ((2.*stry_n + strx_n) * Phi(2*(2*(jphi-1)+iphi),qp,i,j))
            v_diag_qp(iphi,jphi,qp) = v_diag_qp(iphi,jphi,qp) + &
              jac_wt * CS%newton_visc_factor(i,j,qpv) * dstrain_diag_v**2
          endif

          if (grounded_qp) then
            v_diag_qp(iphi,jphi,qp) = v_diag_qp(iphi,jphi,qp) + jac_wt * basal_coef_qp * phi_m_sq
            if (CS%doing_newton) &
              v_diag_qp(iphi,jphi,qp) = v_diag_qp(iphi,jphi,qp) + jac_wt * drag_newt_qp * v_curr_qp**2 * phi_m_sq
          endif
        endif
      enddo ; enddo
    enddo ; enddo

    !element contribution to SW node (node 1, which sees the current element as element 4)
    u_diag_b(I-1,J-1,4) = 0.25*((u_diag_qp(1,1,1)+u_diag_qp(1,1,4))+(u_diag_qp(1,1,2)+u_diag_qp(1,1,3)))
    v_diag_b(I-1,J-1,4) = 0.25*((v_diag_qp(1,1,1)+v_diag_qp(1,1,4))+(v_diag_qp(1,1,2)+v_diag_qp(1,1,3)))

    !element contribution to NW node (node 3, which sees the current element as element 2)
    u_diag_b(I-1,J  ,2) = 0.25*((u_diag_qp(1,2,1)+u_diag_qp(1,2,4))+(u_diag_qp(1,2,2)+u_diag_qp(1,2,3)))
    v_diag_b(I-1,J  ,2) = 0.25*((v_diag_qp(1,2,1)+v_diag_qp(1,2,4))+(v_diag_qp(1,2,2)+v_diag_qp(1,2,3)))

    !element contribution to SE node (node 2, which sees the current element as element 3)
    u_diag_b(I  ,J-1,3) = 0.25*((u_diag_qp(2,1,1)+u_diag_qp(2,1,4))+(u_diag_qp(2,1,2)+u_diag_qp(2,1,3)))
    v_diag_b(I  ,J-1,3) = 0.25*((v_diag_qp(2,1,1)+v_diag_qp(2,1,4))+(v_diag_qp(2,1,2)+v_diag_qp(2,1,3)))

    !element contribution to NE node (node 4, which sees the current element as element 1)
    u_diag_b(I  ,J  ,1) = 0.25*((u_diag_qp(2,2,1)+u_diag_qp(2,2,4))+(u_diag_qp(2,2,2)+u_diag_qp(2,2,3)))
    v_diag_b(I  ,J  ,1) = 0.25*((v_diag_qp(2,2,1)+v_diag_qp(2,2,4))+(v_diag_qp(2,2,2)+v_diag_qp(2,2,3)))

    if (CS%GL_regularize .and. CS%ground_frac(i,j) > 0.0 .and. CS%ground_frac(i,j) < 1.0) then
      ! Subgrid grounding-line: evaluate basal friction diagonal at each grounded sub-quadrature point.
      ! Returns separate u_diag_sub and v_diag_sub (differ in Newton term: u^2 vs v^2).
      ! The sub-qp flotation test handles grounding fraction; no external ground_frac scaling needed.
      Hcell(:,:) = H_node(I-1:I,J-1:J)
      if (do_DG) then
        ! h_nodal_cell is used inside only as a flotation measure (sub-qp gate +
        ! effective pressure), so the gate field is passed (h_flot under
        ! DG_GL_GATE_CONTINUOUS).
        call CG_diagonal_subgrid_basal(CS, G, US, Phisub, Hcell, &
            u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
            CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
            G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
            use_DG=.true., h_shelf_cell=h_shelf(i,j), &
            h_nodal_cell=hgate(i,j,:,:), &
            bed_corners=CS%bed_node(I-1:I,J-1:J))
      else
        call CG_diagonal_subgrid_basal(CS, G, US, Phisub, Hcell, &
            u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
            CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
            G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j))
      endif
      if (CS%umask(I-1,J-1)==1) u_diag_b(I-1,J-1,4) = u_diag_b(I-1,J-1,4) + u_diag_sub(1,1)
      if (CS%umask(I-1,J  )==1) u_diag_b(I-1,J  ,2) = u_diag_b(I-1,J  ,2) + u_diag_sub(1,2)
      if (CS%umask(I  ,J-1)==1) u_diag_b(I  ,J-1,3) = u_diag_b(I  ,J-1,3) + u_diag_sub(2,1)
      if (CS%umask(I  ,J  )==1) u_diag_b(I  ,J  ,1) = u_diag_b(I  ,J  ,1) + u_diag_sub(2,2)
      if (CS%vmask(I-1,J-1)==1) v_diag_b(I-1,J-1,4) = v_diag_b(I-1,J-1,4) + v_diag_sub(1,1)
      if (CS%vmask(I-1,J  )==1) v_diag_b(I-1,J  ,2) = v_diag_b(I-1,J  ,2) + v_diag_sub(1,2)
      if (CS%vmask(I  ,J-1)==1) v_diag_b(I  ,J-1,3) = v_diag_b(I  ,J-1,3) + v_diag_sub(2,1)
      if (CS%vmask(I  ,J  )==1) v_diag_b(I  ,J  ,1) = v_diag_b(I  ,J  ,1) + v_diag_sub(2,2)
    endif
  endif ; enddo ; enddo

  do J=jsc-2,jec+1 ; do I=isc-2,iec+1
    u_diagonal(I,J) = (u_diag_b(I,J,1)+u_diag_b(I,J,4)) + (u_diag_b(I,J,2)+u_diag_b(I,J,3))
    v_diagonal(I,J) = (v_diag_b(I,J,1)+v_diag_b(I,J,4)) + (v_diag_b(I,J,2)+v_diag_b(I,J,3))
  enddo ; enddo

end subroutine matrix_diagonal

!> Compute subgrid grounding-line basal traction contributions for the preconditioner diagonal.
!! Evaluates friction at each grounded sub-quadrature point. Returns separate u and v diagonals
!! because the Newton term uses u^2 for the u-block and v^2 for the v-block.
!! The sub-qp flotation test handles partial grounding; no external ground_frac scaling needed.
subroutine CG_diagonal_subgrid_basal(CS, G, US, Phisub, H_node, U_curr, V_curr, &
                                     bathyT, dens_ratio, i_elem, j_elem, fB_e, u_diag, v_diag, &
                                     dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT, &
                                     use_DG, h_shelf_cell, h_nodal_cell, bed_corners)
  type(ice_shelf_dyn_CS), intent(in) :: CS      !< Ice shelf control structure
  type(ocean_grid_type),  intent(in) :: G       !< The grid structure
  type(unit_scale_type),  intent(in) :: US      !< Unit conversion factors
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub  !< Sub-grid quadrature weights [nondim]
  real, dimension(2,2),   intent(in) :: H_node  !< Ice thickness at element corners [Z ~> m]
  real, dimension(2,2),   intent(in) :: U_curr  !< Frozen u^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_curr  !< Frozen v^k at element corners [L T-1 ~> m s-1]
  real,                   intent(in) :: bathyT  !< Ocean bathymetry depth at tracer point [Z ~> m]
  real,                   intent(in) :: dens_ratio !< Ice density / water density [nondim]
  integer,                intent(in) :: i_elem  !< Tracer-grid i-index of the element
  integer,                intent(in) :: j_elem  !< Tracer-grid j-index of the element
  real,                   intent(in) :: fB_e    !< Element Coulomb parameter fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  real, dimension(2,2),   intent(out) :: u_diag !< Nodal u-diagonal entries [R L2 Z T-1 ~> kg s-1]
  real, dimension(2,2),   intent(out) :: v_diag !< Nodal v-diagonal entries [R L2 Z T-1 ~> kg s-1]
  real,                   intent(in)  :: dxCv_S !< The cell width at the southern (v-point) edge [L ~> m]
  real,                   intent(in)  :: dxCv_N !< The cell width at the northern (v-point) edge [L ~> m]
  real,                   intent(in)  :: dyCu_W !< The cell height at the western (u-point) edge [L ~> m]
  real,                   intent(in)  :: dyCu_E !< The cell height at the eastern (u-point) edge [L ~> m]
  real,                   intent(in)  :: IareaT !< The inverse of the cell area at the tracer point [L-2 ~> m-2]
  logical,       optional, intent(in) :: use_DG       !< If true, use DG thickness and bed_node [nondim]
  real,          optional, intent(in) :: h_shelf_cell  !< Cell-averaged ice thickness [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: h_nodal_cell !< Q1 thickness at the 4 cell corners, used
                                          !! only as a flotation measure (sub-qp gate + effective
                                          !! pressure); the caller passes h_flot here under
                                          !! DG_GL_GATE_CONTINUOUS [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: bed_corners !< Bed elevation at element corners [Z ~> m]

  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: u_diag_sub, v_diag_sub
  real, dimension(2,2,2,2) :: u_diag_qp_nd, v_diag_qp_nd  ! Per-qp nodal diagonal entries (qx,qy,m,n),
                                                         ! pair-summed for rotation invariance
  real :: hloc           ! Local sub-cell ice thickness [Z ~> m]
  real :: bed_sub        ! Bed elevation at sub-quadrature point [Z ~> m]
  real :: u_curr_loc     ! Frozen u^k interpolated to sub-qp [L T-1 ~> m s-1]
  real :: v_curr_loc     ! Frozen v^k interpolated to sub-qp [L T-1 ~> m s-1]
  real :: unorm2_loc     ! Regularized |u^k|^2 at sub-qp [L2 T-2 ~> m2 s-2]
  real :: basal_coef_loc ! Picard friction coefficient at sub-qp [R L2 Z T-1 ~> kg s-1]
  real :: drag_newt_loc  ! Newton drag coefficient at sub-qp [R Z T ~> kg m-2 s]
  real :: phi_mn_sq      ! Squared basis function value at sub-qp [nondim]
  real :: contrib        ! Quadrature weight contribution [nondim]
  real :: coef_prefactor ! Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real :: min_trac_area  ! Minimum area-integrated traction floor [R L2 Z T-1 ~> kg s-1]
  real :: eps_vel2       ! Velocity regularization squared [L2 T-2 ~> m2 s-2]
  real :: jac_sub_wt ! Per-sub-cell-QP metric correction |J_sub|/areaT [nondim]
  real :: a, d      ! Interpolated cell-edge spacings at the sub-cell QP [L ~> m]
  real :: subarea        ! Fractional sub-cell area [nondim]
  real :: fB_local       ! Coulomb fB at sub-qp (DG mode) [(T L-1)^CF_PostPeak]
  real :: rho_oi_ratio   ! density_ocean / density_ice [nondim]
  real :: rho_ice_g_LtoZ ! US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  real :: xi_sub, eta_sub ! DG reference coords at sub-qp ([-0.5,0.5]) [nondim]
  logical :: do_DG       ! Local flag for DG mode
  integer :: nsub, i, j, qx, qy, m, n

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  coef_prefactor = CS%coef_prefactor(i_elem,j_elem)
  min_trac_area  = CS%min_basal_traction * G%areaT(i_elem,j_elem)
  eps_vel2 = CS%eps_glen_min**2 * ((G%dxT(i_elem,j_elem)**2) + (G%dyT(i_elem,j_elem)**2))

  do_DG = .false.
  if (present(use_DG)) do_DG = use_DG
  if (do_DG) then
    rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
    rho_ice_g_LtoZ = US%L_to_Z * CS%density_ice * CS%g_Earth
  endif

  u_diag_sub(:,:,:,:) = 0.0 ; v_diag_sub(:,:,:,:) = 0.0

  do j=1,nsub ; do i=1,nsub
    ! Zero the 4-qp per-node buffer so ungrounded qp contribute exactly 0.
    u_diag_qp_nd(:,:,:,:) = 0.0 ; v_diag_qp_nd(:,:,:,:) = 0.0
    do qy=1,2 ; do qx=1,2
      if (do_DG) then
        ! Nodal Q1 evaluation at the sub-QP via Phisub corner-basis weights.
        hloc = ((Phisub(qx,qy,i,j,1,1)*h_nodal_cell(1,1)) + (Phisub(qx,qy,i,j,2,2)*h_nodal_cell(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*h_nodal_cell(1,2)) + (Phisub(qx,qy,i,j,2,1)*h_nodal_cell(2,1)))
        hloc = max(hloc, CS%min_h_shelf)
        bed_sub = ((Phisub(qx,qy,i,j,1,1)*bed_corners(1,1)) + (Phisub(qx,qy,i,j,2,2)*bed_corners(2,2))) + &
                  ((Phisub(qx,qy,i,j,1,2)*bed_corners(1,2)) + (Phisub(qx,qy,i,j,2,1)*bed_corners(2,1)))
      else
        hloc = ((Phisub(qx,qy,i,j,1,1)*H_node(1,1)) + (Phisub(qx,qy,i,j,2,2)*H_node(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*H_node(1,2)) + (Phisub(qx,qy,i,j,2,1)*H_node(2,1)))
        bed_sub = bathyT
      endif

      if (dens_ratio * hloc - bed_sub > 0) then  ! grounded sub-qp
        u_curr_loc = (((Phisub(qx,qy,i,j,1,1)*U_curr(1,1)) + (Phisub(qx,qy,i,j,2,2)*U_curr(2,2))) + &
                      ((Phisub(qx,qy,i,j,1,2)*U_curr(1,2)) + (Phisub(qx,qy,i,j,2,1)*U_curr(2,1))))
        v_curr_loc = (((Phisub(qx,qy,i,j,1,1)*V_curr(1,1)) + (Phisub(qx,qy,i,j,2,2)*V_curr(2,2))) + &
                      ((Phisub(qx,qy,i,j,1,2)*V_curr(1,2)) + (Phisub(qx,qy,i,j,2,1)*V_curr(2,1))))

        unorm2_loc = ((u_curr_loc**2) + (v_curr_loc**2)) + eps_vel2

        ! Compute Coulomb fB at this sub-qp when in DG mode
        if (do_DG .and. CS%CoulombFriction) then
          fB_local = compute_fB_local(hloc, bed_sub, rho_oi_ratio, rho_ice_g_LtoZ, &
              CS%C_basal_friction(i_elem,j_elem), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
              CS%CF_PostPeak, CS%n_basal_fric)
        else
          fB_local = fB_e
        endif

        call compute_basal_coef(unorm2_loc, coef_prefactor, min_trac_area, fB_local, &
            CS%n_basal_fric, CS%CoulombFriction, CS%CF_PostPeak, US%L_T_to_m_s, .true., &
            basal_coef_loc, drag_newt_loc)
        ! Interpolate cell-edge metrics to the sub-cell QP using the bilinear shape function values
        ! from bilinear_shape_functions_subgrid.  Marginal sums of Phisub give the interpolation
        ! weights: sum over k=1 nodes gives (1-y); k=2 gives y; l=1 gives (1-x); l=2 gives x.
        ! This is analogous to jac_wt = CS%Jac(qp,i,j) * G%IareaT(i,j) in the regular routines.
        a = (dxCv_S * (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1))) + &  ! (1-y) * dxCv_S
            (dxCv_N * (Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)))      !  + y  * dxCv_N
        d = (dyCu_W * (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2))) + &  ! (1-x) * dyCu_W
            (dyCu_E * (Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)))      !  + x  * dyCu_E
        jac_sub_wt = 0.25 * subarea * (a * d) * IareaT

        do n=1,2 ; do m=1,2
          phi_mn_sq = Phisub(qx,qy,i,j,m,n)**2
          contrib   = jac_sub_wt * phi_mn_sq
          if (CS%doing_newton) then
            u_diag_qp_nd(qx,qy,m,n) = contrib * (basal_coef_loc + drag_newt_loc * u_curr_loc**2)
            v_diag_qp_nd(qx,qy,m,n) = contrib * (basal_coef_loc + drag_newt_loc * v_curr_loc**2)
          else
            u_diag_qp_nd(qx,qy,m,n) = contrib * basal_coef_loc
            v_diag_qp_nd(qx,qy,m,n) = contrib * basal_coef_loc
          endif
        enddo ; enddo
      endif
    enddo ; enddo

    do n=1,2 ; do m=1,2
      u_diag_sub(i,j,m,n) = (u_diag_qp_nd(1,1,m,n) + u_diag_qp_nd(2,2,m,n)) + &
                            (u_diag_qp_nd(1,2,m,n) + u_diag_qp_nd(2,1,m,n))
      v_diag_sub(i,j,m,n) = (v_diag_qp_nd(1,1,m,n) + v_diag_qp_nd(2,2,m,n)) + &
                            (v_diag_qp_nd(1,2,m,n) + v_diag_qp_nd(2,1,m,n))
    enddo ; enddo
  enddo ; enddo

  do n=1,2 ; do m=1,2
    call sum_square_matrix(u_diag(m,n), u_diag_sub(:,:,m,n), nsub)
    call sum_square_matrix(v_diag(m,n), v_diag_sub(:,:,m,n), nsub)
  enddo ; enddo

end subroutine CG_diagonal_subgrid_basal

!> Post_data calls related to ice-sheet flux divergence, strain-rate, and deviatoric stress
subroutine IS_dynamics_post_data_2(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJB_(G)) :: H_node ! Ice shelf thickness at corners [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Hu  ! Ice shelf u_flux at corners [Z L T-1 ~> m2 s-1].
  real, dimension(SZDIB_(G),SZDJB_(G)) :: Hv  ! Ice shelf v_flux at corners [Z L T-1 ~> m2 s-1].
  real, dimension(SZDI_(G),SZDJ_(G)) :: Hux  ! Ice shelf d(u_flux)/dx at cell centers [Z T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)) :: Hvy  ! Ice shelf d(v_flux)/dy at cell centers [Z T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)) :: flux_div ! horizontal flux divergence div(uH) [Z T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G),3) :: strain_rate ! strain-rate components xx,yy, and xy [T-1 ~> s-1]
  real, dimension(SZDI_(G),SZDJ_(G),2) :: p_strain_rate ! horizontal principal strain-rates [T-1 ~> s-1]
  real, dimension(SZDI_(G),SZDJ_(G),3) :: dev_stress ! deviatoric stress components xx,yy, and xy [R L Z T-2 ~> Pa]
  real, dimension(SZDI_(G),SZDJ_(G),2) :: p_dev_stress ! horizontal principal deviatoric stress [R L Z T-2 ~> Pa]
  real, dimension(SZDI_(G),SZDJ_(G))  :: ice_visc ! area-averaged ice viscosity [R L2 T-1 ~> Pa s]
  real :: p1,p2 ! Used to calculate strain-rate principal components [T-1 ~> s-1]
  integer :: i, j

  !Allocate the gradient basis functions for 1 cell-centered quadrature point per cell
  if (.not. associated(CS%PhiC)) then
    allocate(CS%PhiC(1:8,G%isc:G%iec,G%jsc:G%jec), source=0.0)
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      call bilinear_shape_fn_grid_1qp(G, i, j, CS%PhiC(:,i,j))
    enddo ; enddo
  endif

  !Calculate flux divergence and its components
  if (CS%id_duHdx > 0 .or. CS%id_dvHdy > 0 .or. CS%id_fluxdiv > 0) then
    if (CS%use_DG_thickness) then
      call interpolate_H_to_B_DG(G, ISS%h_shelf, CS%h_nodal, ISS%hmask, &
                                 H_node, CS%min_h_shelf)
    else
      call interpolate_H_to_B(G, ISS%h_shelf, ISS%hmask, H_node, CS%min_h_shelf)
    endif

    Hu(:,:) = 0.0 ; Hv(:,:) = 0.0 ; Hux(:,:) = 0.0 ; Hvy(:,:) = 0.0 ; flux_div(:,:) = 0.0
    do J=G%jscB,G%jecB ; do I=G%iscB,G%iecB
      if (CS%umask(I,J) > 0) then
        Hu(I,J) = (H_node(I,J) * CS%u_shelf(I,J))
      endif
      if (CS%vmask(I,J) > 0) then
        Hv(I,J) = (H_node(I,J) * CS%v_shelf(I,J))
      endif
    enddo ; enddo

    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 3)) then
        !components of flux divergence at cell centers
        Hux(i,j) = (((Hu(I-1,J-1) * CS%PhiC(1,i,j)) + (Hu(I,J  ) * CS%PhiC(7,i,j))) + &
                    ((Hu(I-1,J  ) * CS%PhiC(5,i,j)) + (Hu(I,J-1) * CS%PhiC(3,i,j))))

        Hvy(i,j) = (((Hv(I-1,J-1) * CS%PhiC(2,i,j)) + (Hv(I,J  ) * CS%PhiC(8,i,j))) + &
                    ((Hv(I-1,J  ) * CS%PhiC(6,i,j)) + (Hv(I,J-1) * CS%PhiC(4,i,j))))
        flux_div(i,j) = Hux(i,j) + Hvy(i,j)
      endif
    enddo ; enddo

    if (CS%id_duHdx > 0)   call post_data(CS%id_duHdx, Hux, CS%diag)
    if (CS%id_dvHdy > 0)   call post_data(CS%id_dvHdy, Hvy, CS%diag)
    if (CS%id_fluxdiv > 0) call post_data(CS%id_fluxdiv, flux_div, CS%diag)
  endif

  if (CS%id_devstress_xx > 0  .or. CS%id_devstress_yy > 0  .or. CS%id_devstress_xy > 0  .or. &
      CS%id_strainrate_xx > 0 .or. CS%id_strainrate_yy > 0 .or. CS%id_strainrate_xy > 0 .or. &
      CS%id_pdevstress_1 > 0  .or. CS%id_pdevstress_2 > 0  .or. &
      CS%id_pstrainrate_1 > 0 .or. CS%id_pstrainrate_2 > 0) then

    strain_rate(:,:,:) = 0.0
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      !strain-rates at cell centers
      if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 3)) then
        !strain_rate(:,:,1) = strain_rate_xx(:,:) = ux(:,:)
        strain_rate(i,j,1) = (((CS%u_shelf(I-1,J-1) * CS%PhiC(1,i,j)) + (CS%u_shelf(I,J  ) * CS%PhiC(7,i,j))) + &
                              ((CS%u_shelf(I-1,J  ) * CS%PhiC(5,i,j)) + (CS%u_shelf(I,J-1) * CS%PhiC(3,i,j))))
        !strain_rate(:,:,2) = strain_rate_yy(:,:) = uy(:,:)
        strain_rate(i,j,2) = (((CS%v_shelf(I-1,J-1) * CS%PhiC(2,i,j)) + (CS%v_shelf(I,J  ) * CS%PhiC(8,i,j))) + &
                              ((CS%v_shelf(I-1,J  ) * CS%PhiC(6,i,j)) + (CS%v_shelf(I,J-1) * CS%PhiC(4,i,j))))
        !strain_rate(:,:,3) = strain_rate_xy(:,:) = 0.5 * (uy(:,:) + vy(:,:))
        strain_rate(i,j,3) = 0.5 * ((((CS%u_shelf(I-1,J-1) * CS%PhiC(2,i,j)) + (CS%u_shelf(I,J  ) * CS%PhiC(8,i,j))) + &
                                     ((CS%u_shelf(I-1,J  ) * CS%PhiC(6,i,j)) + (CS%u_shelf(I,J-1) * CS%PhiC(4,i,j))))+ &
                                    (((CS%v_shelf(I-1,J-1) * CS%PhiC(1,i,j)) + (CS%v_shelf(I,J  ) * CS%PhiC(7,i,j))) + &
                                     ((CS%v_shelf(I-1,J  ) * CS%PhiC(5,i,j)) + (CS%v_shelf(I,J-1) * CS%PhiC(3,i,j)))))
      endif
    enddo ; enddo


    if (CS%id_strainrate_xx > 0) call post_data(CS%id_strainrate_xx, strain_rate(:,:,1), CS%diag)
    if (CS%id_strainrate_yy > 0) call post_data(CS%id_strainrate_yy, strain_rate(:,:,2), CS%diag)
    if (CS%id_strainrate_xy > 0) call post_data(CS%id_strainrate_xy, strain_rate(:,:,3), CS%diag)

    if (CS%id_pstrainrate_1 > 0 .or. CS%id_pstrainrate_2 > 0 .or. &
        CS%id_pdevstress_1  > 0 .or. CS%id_pdevstress_2  > 0) then
      p_strain_rate(:,:,:) = 0.0
      do j=G%jsc,G%jec ; do i=G%isc,G%iec
        p1 = 0.5*( strain_rate(i,j,1) + strain_rate(i,j,2))
        p2 = sqrt( (( 0.5 * (strain_rate(i,j,1) - strain_rate(i,j,2)) )**2) + (strain_rate(i,j,3)**2) )
        p_strain_rate(i,j,1) = p1+p2 !Max horizontal principal strain-rate
        p_strain_rate(i,j,2) = p1-p2 !Min horizontal principal strain-rate
      enddo ; enddo

      if (CS%id_pstrainrate_1 > 0) call post_data(CS%id_pstrainrate_1, p_strain_rate(:,:,1), CS%diag)
      if (CS%id_pstrainrate_2 > 0) call post_data(CS%id_pstrainrate_2, p_strain_rate(:,:,2), CS%diag)
    endif

    if (CS%id_devstress_xx > 0 .or. CS%id_devstress_yy > 0 .or. CS%id_devstress_xy > 0 .or. &
        CS%id_pdevstress_1 > 0 .or. CS%id_pdevstress_2 > 0) then

      call ice_visc_diag(CS,G,ice_visc)

      if (CS%id_devstress_xx > 0 .or. CS%id_devstress_yy > 0 .or. CS%id_devstress_xy > 0) then
        dev_stress(:,:,:)=0.0
        do j=G%jsc,G%jec ; do i=G%isc,G%iec
          if (ISS%h_shelf(i,j)>0) then
            dev_stress(i,j,1) = 2*ice_visc(i,j)*strain_rate(i,j,1)/ISS%h_shelf(i,j) !deviatoric stress xx
            dev_stress(i,j,2) = 2*ice_visc(i,j)*strain_rate(i,j,2)/ISS%h_shelf(i,j) !deviatoric stress yy
            dev_stress(i,j,3) = 2*ice_visc(i,j)*strain_rate(i,j,3)/ISS%h_shelf(i,j) !deviatoric stress xy
          endif
        enddo ; enddo
        if (CS%id_devstress_xx > 0) call post_data(CS%id_devstress_xx, dev_stress(:,:,1), CS%diag)
        if (CS%id_devstress_yy > 0) call post_data(CS%id_devstress_yy, dev_stress(:,:,2), CS%diag)
        if (CS%id_devstress_xy > 0) call post_data(CS%id_devstress_xy, dev_stress(:,:,3), CS%diag)
      endif

      if (CS%id_pdevstress_1 > 0 .or. CS%id_pdevstress_2 > 0) then
        p_dev_stress(:,:,:)=0.0
        do j=G%jsc,G%jec ; do i=G%isc,G%iec
          if (ISS%h_shelf(i,j)>0) then
            p_dev_stress(i,j,1) = 2*ice_visc(i,j)*p_strain_rate(i,j,1)/ISS%h_shelf(i,j) !max horiz principal dev stress
            p_dev_stress(i,j,2) = 2*ice_visc(i,j)*p_strain_rate(i,j,2)/ISS%h_shelf(i,j) !min horiz principal dev stress
          endif
        enddo ; enddo
        if (CS%id_pdevstress_1 > 0) call post_data(CS%id_pdevstress_1, p_dev_stress(:,:,1), CS%diag)
        if (CS%id_pdevstress_2 > 0) call post_data(CS%id_pdevstress_2, p_dev_stress(:,:,2), CS%diag)
      endif
    endif
  endif
end subroutine IS_dynamics_post_data_2

!> Update depth integrated viscosity, based on horizontal strain rates
subroutine calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: u_shlf !< The zonal ice shelf velocity [L T-1 ~> m s-1].
  real, dimension(G%IsdB:G%IedB,G%JsdB:G%JedB), &
                          intent(inout) :: v_shlf !< The meridional ice shelf velocity [L T-1 ~> m s-1].

! update DEPTH_INTEGRATED viscosity, based on horizontal strain rates - this is for bilinear FEM solve


! this may be subject to change later... to make it "hybrid"
!  real, dimension(SZDIB_(G),SZDJB_(G)) ::  eII, ux, uy, vx, vy
  integer :: i, j, iscq, iecq, jscq, jecq, isd, jsd, ied, jed, iegq, jegq, iq, jq
  integer :: giec, gjec, gisc, gjsc, isc, jsc, iec, jec, is, js
  real :: Visc_coef, n_g
  real :: ux, uy, vx, vy
  real :: eps_min   ! Velocity shears [T-1 ~> s-1]
  real :: h_gp      ! DG-evaluated ice thickness at quadrature point [Z ~> m]
  real, dimension(2) :: xquad ! Gauss quadrature positions on [0,1] [nondim]
  logical :: model_qp1, model_qp4

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%domain%nihalo+1 ; gjsc = G%domain%njhalo+1
  giec = G%domain%niglobal+gisc ; gjec = G%domain%njglobal+gjsc
  is = iscq - 1 ; js = jscq - 1

  if (trim(CS%ice_viscosity_compute) == "MODEL") then
    if (CS%visc_qps==1) then
      model_qp1=.true.
      model_qp4=.false.
    else
      model_qp1=.false.
      model_qp4=.true.
    endif
  endif

  n_g = CS%n_glen ; eps_min = CS%eps_glen_min
  xquad(1) = 0.5 * (1.0 - sqrt(1.0/3.0)) ; xquad(2) = 0.5 * (1.0 + sqrt(1.0/3.0))

  do j=jsc,jec ; do i=isc,iec

    if ((ISS%hmask(i,j) == 1) .OR. (ISS%hmask(i,j) == 3)) then

      if (trim(CS%ice_viscosity_compute) == "CONSTANT") then
        CS%ice_visc(i,j,1) = 1e15 * (US%kg_m3_to_R*US%m_to_L*US%m_s_to_L_T) * &
                             (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf))
        ! constant viscocity for debugging
      elseif (trim(CS%ice_viscosity_compute) == "OBS") then
        if (CS%AGlen_visc(i,j) >0) then
          CS%ice_visc(i,j,1) = (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf)) * &
                               max(CS%AGlen_visc(i,j) ,CS%min_ice_visc)
        endif
        ! Here CS%Aglen_visc(i,j) is the ice viscosity [R L2 T-1 ~> Pa s] computed from obs and read from a file
      elseif (model_qp1) then
        ! calculate viscosity at 1 cell-centered quadrature point per cell

        Visc_coef = (CS%AGlen_visc(i,j))**(-1./n_g)
        ! Units of Aglen_visc [Pa-(n_g) s-1]

        ux = ((u_shlf(I-1,J-1) * CS%PhiC(1,i,j)) + &
              (u_shlf(I,J) * CS%PhiC(7,i,j))) + &
             ((u_shlf(I-1,J) * CS%PhiC(5,i,j)) + &
              (u_shlf(I,J-1) * CS%PhiC(3,i,j)))

        vx = ((v_shlf(I-1,J-1) * CS%PhiC(1,i,j)) + &
              (v_shlf(I,J) * CS%PhiC(7,i,j))) + &
             ((v_shlf(I-1,J) * CS%PhiC(5,i,j)) + &
              (v_shlf(I,J-1) * CS%PhiC(3,i,j)))

        uy = ((u_shlf(I-1,J-1) * CS%PhiC(2,i,j)) + &
              (u_shlf(I,J) * CS%PhiC(8,i,j))) + &
             ((u_shlf(I-1,J) * CS%PhiC(6,i,j)) + &
              (u_shlf(I,J-1) * CS%PhiC(4,i,j)))

        vy = ((v_shlf(I-1,J-1) * CS%PhiC(2,i,j)) + &
              (v_shlf(I,J) * CS%PhiC(8,i,j))) + &
             ((v_shlf(I-1,J) * CS%PhiC(6,i,j)) + &
              (v_shlf(I,J-1) * CS%PhiC(4,i,j)))

        CS%ice_visc(i,j,1) = (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf)) * &
            max(0.5 * Visc_coef * &
            (US%s_to_T**2 * (((ux**2) + (vy**2)) + ((ux*vy) + 0.25*((uy+vx)**2)) + eps_min**2))**((1.-n_g)/(2.*n_g)) * &
            (US%Pa_to_RL2_T2*US%s_to_T),CS%min_ice_visc)  ! Rescale after the fractional power law.
        ! Store Newton tangent stiffness data: strain rates and coefficient for Newton iterations.
        ! The Newton correction coefficient is (1/n-1)/2 * ice_visc / eps_e2,
        ! where eps_e2 = ux^2 + vy^2 + ux*vy + (uy+vx)^2/4 + eps_min^2 [T-2].
        ! It is zero where ice_visc is limited by min_ice_visc (viscosity is not smooth there).
        CS%newton_str_ux(i,j,1) = ux ; CS%newton_str_vy(i,j,1) = vy
        CS%newton_str_sh(i,j,1) = uy + vx
        CS%newton_visc_factor(i,j,1) = 0.0
        if (CS%ice_visc(i,j,1) > CS%min_ice_visc * (G%areaT(i,j) * max(ISS%h_shelf(i,j),CS%min_h_shelf))) then
          CS%newton_visc_factor(i,j,1) = ((1./n_g - 1.) / &
              (((ux**2) + (vy**2)) + ((ux*vy) + 0.25*((uy+vx)**2)) + eps_min**2)) * &
              CS%ice_visc(i,j,1)
        endif
      elseif (model_qp4) then
        !calculate viscosity at 4 quadrature points per cell

        Visc_coef = (CS%AGlen_visc(i,j))**(-1./n_g)

        do iq=1,2 ; do jq=1,2

          ! Evaluate ice thickness at this Gauss point
          if (CS%use_DG_thickness) then
            h_gp = ((CS%h_nodal(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                    (CS%h_nodal(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
                   ((CS%h_nodal(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                    (CS%h_nodal(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
            h_gp = max(h_gp, CS%min_h_shelf)
          else
            h_gp = max(ISS%h_shelf(i,j), CS%min_h_shelf)
          endif

          ux = ((u_shlf(I-1,J-1) * CS%Phi(1,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I,J) * CS%Phi(7,2*(jq-1)+iq,i,j))) + &
               ((u_shlf(I,J-1) * CS%Phi(3,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I-1,J) * CS%Phi(5,2*(jq-1)+iq,i,j)))

          vx = ((v_shlf(I-1,J-1) * CS%Phi(1,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I,J) * CS%Phi(7,2*(jq-1)+iq,i,j))) + &
               ((v_shlf(I,J-1) * CS%Phi(3,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I-1,J) * CS%Phi(5,2*(jq-1)+iq,i,j)))

          uy = ((u_shlf(I-1,J-1) * CS%Phi(2,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I,J) * CS%Phi(8,2*(jq-1)+iq,i,j))) + &
               ((u_shlf(I,J-1) * CS%Phi(4,2*(jq-1)+iq,i,j)) + &
                (u_shlf(I-1,J) * CS%Phi(6,2*(jq-1)+iq,i,j)))

          vy = ((v_shlf(I-1,J-1) * CS%Phi(2,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I,J) * CS%Phi(8,2*(jq-1)+iq,i,j))) + &
               ((v_shlf(I,J-1) * CS%Phi(4,2*(jq-1)+iq,i,j)) + &
                (v_shlf(I-1,J) * CS%Phi(6,2*(jq-1)+iq,i,j)))

          CS%ice_visc(i,j,2*(jq-1)+iq) = (G%areaT(i,j) * h_gp) * &
              max(0.5 * Visc_coef * &
              (US%s_to_T**2*(((ux**2) + (vy**2)) + ((ux*vy) + 0.25*((uy+vx)**2)) + eps_min**2))**((1.-n_g)/(2.*n_g)) * &
              (US%Pa_to_RL2_T2*US%s_to_T),CS%min_ice_visc)  ! Rescale after the fractional power law.
          ! Store Newton tangent stiffness data at each quadrature point.
          CS%newton_str_ux(i,j,2*(jq-1)+iq) = ux ; CS%newton_str_vy(i,j,2*(jq-1)+iq) = vy
          CS%newton_str_sh(i,j,2*(jq-1)+iq) = (uy + vx)
          CS%newton_visc_factor(i,j,2*(jq-1)+iq) = 0.0
          if (CS%ice_visc(i,j,2*(jq-1)+iq) > CS%min_ice_visc * (G%areaT(i,j) * h_gp)) then
            CS%newton_visc_factor(i,j,2*(jq-1)+iq) = ((1./n_g - 1.) / &
                (((ux**2) + (vy**2)) + ((ux*vy) + 0.25*((uy+vx)**2)) + eps_min**2)) * &
                CS%ice_visc(i,j,2*(jq-1)+iq)
          endif
        enddo ; enddo
      endif
    endif
  enddo ; enddo

end subroutine calc_shelf_visc

!> Pre-compute element-level basal friction prefactors for quadrature-point evaluation.
subroutine calc_shelf_basal_prefactors(CS, ISS, G, US)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< Ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< Ice shelf state (hmask, h_shelf)
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure
  type(unit_scale_type),  intent(in)    :: US  !< Unit conversion factors

  integer :: i, j, isd, ied, jsd, jed
  real :: Hf  ! Floatation thickness [Z ~> m]
  real :: fN  ! Effective pressure for Coulomb friction [R Z L T-2 ~> Pa]

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j = jsd, jed ; do i = isd, ied
    CS%coef_prefactor(i,j) = G%areaT(i,j) * CS%C_basal_friction(i,j) * US%L_T_to_m_s
    if (CS%use_DG_thickness) then
      ! DG mode: fB is computed at each quadrature point via compute_fB_local; zero the element value.
      CS%fB_elem(i,j) = 0.0
    elseif (CS%CoulombFriction .and. (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3)) then
      Hf = max((CS%density_ocean_avg/CS%density_ice) * CS%bed_elev(i,j), 0.0)
      fN = max((US%L_to_Z*(CS%density_ice * CS%g_Earth) * &
                (max(ISS%h_shelf(i,j), CS%min_h_shelf) - Hf)), CS%CF_MinN)
      CS%fB_elem(i,j) = CS%alpha_coulomb * &
          (CS%C_basal_friction(i,j) / (CS%CF_Max * fN))**(CS%CF_PostPeak/CS%n_basal_fric)
    else
      CS%fB_elem(i,j) = 0.0
    endif
  enddo ; enddo

end subroutine calc_shelf_basal_prefactors

!> Compute area-averaged basal shear stress [R L T-1 ~> Pa s m-1] and return it in basal_tr.
!! Uses CS%u_shelf and CS%v_shelf for velocities and G%US for unit conversions.
subroutine calc_shelf_taub(CS, ISS, G, basal_tr)
  type(ice_shelf_dyn_CS), intent(in)  :: CS  !< Ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)  :: ISS !< A structure with elements that describe
                                             !! the ice-shelf state
  type(ocean_grid_type),  intent(in)  :: G   !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(out) :: basal_tr !< Area-averaged basal traction [R L T-1 ~> Pa s m-1]

  integer :: i, j
  real :: umid, vmid    ! Cell-center velocity averages [L T-1 ~> m s-1]
  real :: eps_min       ! Minimal strain rate [T-1 ~> s-1]
  real :: unorm         ! Velocity magnitude in mks units [m s-1]
  real :: alpha         ! Coulomb coefficient [nondim]
  real :: Hf            ! Floatation thickness for Coulomb friction [Z ~> m]
  real :: fN            ! Effective pressure for Coulomb friction [R Z L T-2 ~> Pa]
  real :: fB            ! Coulomb friction factor [(T L-1)^CS%CF_PostPeak]
  real :: fBuq          ! fB * unorm^CF_PostPeak [nondim]
  real :: unorm_code2   ! Squared velocity magnitude in code units [L2 T-2 ~> m2 s-2]
  real :: basal_trac    ! Area-integrated traction coefficient [R Z L2 T-1 ~> kg s-1]

  eps_min = CS%eps_glen_min

  if (CS%CoulombFriction) then
    if (CS%CF_PostPeak /= 1.0) then
      alpha = (CS%CF_PostPeak-1.0)**(CS%CF_PostPeak-1.0) / CS%CF_PostPeak**CS%CF_PostPeak
    else
      alpha = 1.0
    endif
  endif

  basal_tr(:,:) = 0.0

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((ISS%hmask(i,j) == 1) .OR. (ISS%hmask(i,j) == 3)) then
      umid = ((CS%u_shelf(I,J) + CS%u_shelf(I-1,J-1)) + (CS%u_shelf(I,J-1) + CS%u_shelf(I-1,J))) * 0.25
      vmid = ((CS%v_shelf(I,J) + CS%v_shelf(I-1,J-1)) + (CS%v_shelf(I,J-1) + CS%v_shelf(I-1,J))) * 0.25
      unorm_code2 = ((umid**2) + (vmid**2)) + (eps_min**2 * ((G%dxT(i,j)**2) + (G%dyT(i,j)**2)))
      unorm = G%US%L_T_to_m_s * sqrt(unorm_code2)

      !Coulomb friction (Schoof 2005, Gagliardini et al 2007)
      if (CS%CoulombFriction) then
        !Effective pressure
        Hf = max((CS%density_ocean_avg/CS%density_ice) * CS%bed_elev(i,j), 0.0)
        fN = max((G%US%L_to_Z*(CS%density_ice * CS%g_Earth) * (max(ISS%h_shelf(i,j),CS%min_h_shelf) - Hf)), CS%CF_MinN)
        fB = alpha * (CS%C_basal_friction(i,j) / (CS%CF_Max * fN))**(CS%CF_PostPeak/CS%n_basal_fric)
        fBuq = fB * unorm**CS%CF_PostPeak
        basal_trac = ((G%areaT(i,j) * CS%C_basal_friction(i,j)) * &
            (unorm**(CS%n_basal_fric-1.0) / (1.0 + fBuq)**(CS%n_basal_fric))) * &
            G%US%L_T_to_m_s   ! Restore the scaling after the fractional power law.
      else
        !linear (CS%n_basal_fric = 1) or "Weertman"/power-law (CS%n_basal_fric /= 1)
        basal_trac = ((G%areaT(i,j) * CS%C_basal_friction(i,j)) * (unorm**(CS%n_basal_fric-1))) * &
                     G%US%L_T_to_m_s ! Rescale after the fractional power law.
      endif

      basal_trac = max(basal_trac, CS%min_basal_traction * G%areaT(i,j))
      basal_tr(i,j) = basal_trac * G%IareaT(i,j) * CS%ground_frac(i,j)
    endif
  enddo ; enddo

end subroutine calc_shelf_taub

subroutine update_OD_ffrac(CS, G, US, ocean_mass, find_avg)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)     :: US !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: ocean_mass !< The mass per unit area of the ocean [R Z ~> kg m-2].
  logical,                intent(in)    :: find_avg !< If true, find the average of OD and ffrac, and
                                              !! reset the underlying running sums to 0.

  integer :: isc, iec, jsc, jec, i, j
  real    :: I_rho_ocean ! A typical specific volume of the ocean [R-1 ~> m3 kg-1]
  real    :: I_counter

  I_rho_ocean = 1.0 / CS%density_ocean_avg

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  do j=jsc,jec ; do i=isc,iec
    CS%OD_rt(i,j) = CS%OD_rt(i,j) + ocean_mass(i,j)*I_rho_ocean
    if (ocean_mass(i,j)*I_rho_ocean > CS%thresh_float_col_depth) then
      CS%ground_frac_rt(i,j) = CS%ground_frac_rt(i,j) + 1.0
    endif
  enddo ; enddo
  CS%OD_rt_counter = CS%OD_rt_counter + 1

  if (find_avg) then
    I_counter = 1.0 / real(CS%OD_rt_counter)
    do j=jsc,jec ; do i=isc,iec
      CS%ground_frac(i,j) = 1.0 - (CS%ground_frac_rt(i,j) * I_counter)
      CS%OD_av(i,j) = CS%OD_rt(i,j) * I_counter

      CS%OD_rt(i,j) = 0.0 ; CS%ground_frac_rt(i,j) = 0.0 ; CS%OD_rt_counter = 0
    enddo ; enddo

    call pass_var(CS%ground_frac, G%domain, complete=.false.)
    call pass_var(CS%OD_av, G%domain, complete=.true.)
  endif

end subroutine update_OD_ffrac

subroutine update_OD_ffrac_uncoupled(CS, G, h_shelf)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf !< the thickness of the ice shelf [Z ~> m].

  integer :: i, j, isd, ied, jsd, jed
  real    :: rhoi_rhow, OD

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j=jsd,jed
    do i=isd,ied
      OD = CS%bed_elev(i,j) - rhoi_rhow * max(h_shelf(i,j),CS%min_h_shelf)
      if (OD >= 0) then
    ! ice thickness does not take up whole ocean column -> floating
        CS%OD_av(i,j) = OD
        CS%ground_frac(i,j) = 0.
      else
        CS%OD_av(i,j) = 0.
        CS%ground_frac(i,j) = 1.
      endif
    enddo
  enddo

end subroutine update_OD_ffrac_uncoupled

!> Build the continuous (C0) flotation-gate thickness field CS%h_flot from the DG
!! nodal thickness. Each B-grid node value is an average of the h_nodal corner
!! values of all adjacent ice cells (hmask 1 or 3), written back into every
!! participating cell's corner slot so that cells sharing a node hold identical
!! values. With DG_GL_GATE_DEFICIT_SCALE <= 0 the average is arithmetic; with a
!! positive scale s each corner is weighted by 1/(|d_k| + s), where
!! d_k = (rho_i/rho_w)*h_k - bed is that corner's flotation deficit, so the side
!! nearer flotation dominates and a large one-sided jump at the grounding-line
!! face cannot drag the gate far into the lighter cell. s -> infinity recovers
!! the arithmetic mean; s -> 0 pins straddling faces exactly at flotation, which
!! reintroduces a dead band (gate insensitive to thickness changes while the
!! face straddles) and should be avoided. The gate is a locator, not a mass
!! field: it commits the flotation crossing to a single position inside any
!! inter-cell thickness jump, so the grounded/floating gates vary continuously
!! as the grounding line migrates through faces. Must be called by all PEs
!! (contains halo exchanges).
subroutine compute_h_flot(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure used by the ice shelf.

  real    :: h_corner(4) ! h_nodal corner values gathered at the node [Z ~> m]
  real    :: w_sum   ! Sum of corner weights [Z-1 ~> m-1] (deficit weighting) or [nondim]
  real    :: hw_sum  ! Weight-thickness sum [nondim] (deficit weighting) or [Z ~> m]
  real    :: w_k     ! Weight of one corner [Z-1 ~> m-1] (deficit weighting) or [nondim]
  real    :: h_avg   ! Node-averaged gate thickness [Z ~> m]
  real    :: bed_n   ! Bed depth at the node [Z ~> m]
  real    :: rhoi_rhow ! Ice/ocean density ratio [nondim]
  real    :: s_def   ! Deficit-weighting regularization scale [Z ~> m]
  integer :: n_cells ! Number of ice cells (hmask 1 or 3) sharing the node
  integer :: i, j, k, isd, ied, jsd, jed
  logical :: vSW, vSE, vNW, vNE ! True if the cell on that side of the node is ice
  logical :: deficit_weighting  ! True if DG_GL_GATE_DEFICIT_SCALE > 0

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  s_def = CS%dg_gl_gate_deficit_scale
  deficit_weighting = (s_def > 0.0)

  ! Ensure h_nodal halo corners are current before averaging across PE edges.
  call pass_corner_field(CS%h_nodal, G)

  CS%h_flot(:,:,:,:) = 0.0

  ! Loop over B-grid nodes (I,J); the 4 surrounding cells are (i,j) to the SW
  ! [corner (2,2)], (i+1,j) to the SE [corner (1,2)], (i,j+1) to the NW
  ! [corner (2,1)], and (i+1,j+1) to the NE [corner (1,1)].
  do j=jsd,jed-1 ; do i=isd,ied-1
    vSW = (ISS%hmask(i  ,j  ) == 1.0 .or. ISS%hmask(i  ,j  ) == 3.0)
    vSE = (ISS%hmask(i+1,j  ) == 1.0 .or. ISS%hmask(i+1,j  ) == 3.0)
    vNW = (ISS%hmask(i  ,j+1) == 1.0 .or. ISS%hmask(i  ,j+1) == 3.0)
    vNE = (ISS%hmask(i+1,j+1) == 1.0 .or. ISS%hmask(i+1,j+1) == 3.0)

    n_cells = 0
    if (vSW) then ; n_cells = n_cells + 1 ; h_corner(n_cells) = CS%h_nodal(i  ,j  ,2,2) ; endif
    if (vSE) then ; n_cells = n_cells + 1 ; h_corner(n_cells) = CS%h_nodal(i+1,j  ,1,2) ; endif
    if (vNW) then ; n_cells = n_cells + 1 ; h_corner(n_cells) = CS%h_nodal(i  ,j+1,2,1) ; endif
    if (vNE) then ; n_cells = n_cells + 1 ; h_corner(n_cells) = CS%h_nodal(i+1,j+1,1,1) ; endif
    if (n_cells == 0) cycle

    if (deficit_weighting) then
      ! Inverse-deficit weighting: corners nearer flotation dominate. The node
      ! (I,J) = (i,j) is the NE corner of cell (i,j), so the bed there is
      ! bed_node(i,j) (single-valued; the bed field is continuous).
      bed_n = CS%bed_node(i,j)
      w_sum = 0.0 ; hw_sum = 0.0
      do k=1,n_cells
        w_k = 1.0 / (abs((rhoi_rhow * h_corner(k)) - bed_n) + s_def)
        w_sum = w_sum + w_k
        hw_sum = hw_sum + (w_k * h_corner(k))
      enddo
      h_avg = hw_sum / w_sum
    else
      hw_sum = 0.0
      do k=1,n_cells
        hw_sum = hw_sum + h_corner(k)
      enddo
      h_avg = hw_sum / real(n_cells)
    endif

    if (vSW) CS%h_flot(i  ,j  ,2,2) = h_avg
    if (vSE) CS%h_flot(i+1,j  ,1,2) = h_avg
    if (vNW) CS%h_flot(i  ,j+1,2,1) = h_avg
    if (vNE) CS%h_flot(i+1,j+1,1,1) = h_avg
  enddo ; enddo

  ! Fill halo corners that this PE's node loop could not reach.
  call pass_corner_field(CS%h_flot, G)

end subroutine compute_h_flot

!> Set CS%ground_frac to the fraction of sub-grid integration points that are
!! grounded. Only active under GL_regularize=True; in that path the sub-grid uses
!! n_sub_regularize x n_sub_regularize sub-cells with 2x2 Gauss points each (the same
!! quadrature stencil that calc_shelf_driving_stress_DG uses for volume integrals in
!! GL-regularize cells). In other cells the input ground_frac is left untouched, so
!! the binary uncoupled value (set by update_OD_ffrac_uncoupled and reinforced inside
!! ice_shelf_solve_outer) or the smooth running-mean value (set by update_OD_ffrac
!! under GL_couple=True) is preserved.
subroutine compute_ground_frac(CS, ISS, G, H_node)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(in)    :: H_node !< Ice shelf thickness at B-grid corners
                                                  !! (set by interpolate_H_to_B; used only
                                                  !! in the non-DG branch) [Z ~> m].

  real :: rhoi_rhow      ! Ice/ocean density ratio [nondim]
  real :: h_ip, bed_ip   ! Ice thickness and bed elevation at a sub-IP [Z ~> m]
  real :: bed_corners(2,2) ! Bed elevation at the 4 B-grid corners of a cell [Z ~> m]
  real :: H_corners(2,2)   ! Ice thickness at the 4 B-grid corners (non-DG path) [Z ~> m]
  real, dimension(:,:,:,:), pointer :: hgate ! Thickness field for the flotation
                         ! test: h_flot under DG_GL_GATE_CONTINUOUS, else h_nodal [Z ~> m]
  integer :: i, j, isub, jsub, iq, jq, n_total, n_grounded
  integer :: isc, iec, jsc, jec

  if (.not. CS%GL_regularize) return

  if (CS%use_DG_thickness) then
    hgate => CS%h_nodal
    if (CS%dg_gl_gate_continuous) hgate => CS%h_flot
  endif

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  n_total = CS%n_sub_regularize * CS%n_sub_regularize * 4
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  do j=jsc,jec ; do i=isc,iec
    if (ISS%hmask(i,j) /= 1 .and. ISS%hmask(i,j) /= 3) cycle

    ! Gather corner values for bed (and H if non-DG) at this cell.
    if (CS%use_DG_thickness) then
      bed_corners(1,1) = CS%bed_node(i-1,j-1) ; bed_corners(2,1) = CS%bed_node(i,j-1)
      bed_corners(1,2) = CS%bed_node(i-1,j  ) ; bed_corners(2,2) = CS%bed_node(i,j  )
    else
      ! Non-DG: bed is cell-constant in the existing GL detection logic; mirror that
      ! by setting all 4 corner values to bed_elev(i,j).
      bed_corners(:,:) = CS%bed_elev(i,j)
      H_corners(1,1) = H_node(i-1,j-1) ; H_corners(2,1) = H_node(i,j-1)
      H_corners(1,2) = H_node(i-1,j  ) ; H_corners(2,2) = H_node(i,j  )
    endif

    n_grounded = 0
    do jsub=1,CS%n_sub_regularize ; do isub=1,CS%n_sub_regularize
      do jq=1,2 ; do iq=1,2
        if (CS%use_DG_thickness) then
          h_ip = ((CS%Phisub(iq,jq,isub,jsub,1,1)*hgate(i,j,1,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,2,2)*hgate(i,j,2,2))) + &
                 ((CS%Phisub(iq,jq,isub,jsub,2,1)*hgate(i,j,2,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,1,2)*hgate(i,j,1,2)))
          h_ip = max(h_ip, CS%min_h_shelf)
        else
          h_ip = ((CS%Phisub(iq,jq,isub,jsub,1,1)*H_corners(1,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,2,2)*H_corners(2,2))) + &
                 ((CS%Phisub(iq,jq,isub,jsub,2,1)*H_corners(2,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,1,2)*H_corners(1,2)))
          h_ip = max(h_ip, CS%min_h_shelf)
        endif

        bed_ip = ((CS%Phisub(iq,jq,isub,jsub,1,1)*bed_corners(1,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,2,2)*bed_corners(2,2))) + &
                 ((CS%Phisub(iq,jq,isub,jsub,2,1)*bed_corners(2,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,1,2)*bed_corners(1,2)))

        if (rhoi_rhow * h_ip - bed_ip > 0.0) n_grounded = n_grounded + 1
      enddo ; enddo
    enddo ; enddo

    CS%ground_frac(i,j) = real(n_grounded) / real(n_total)
  enddo ; enddo

  call pass_var(CS%ground_frac, G%Domain, complete=.true.)

end subroutine compute_ground_frac

subroutine change_in_draft(CS, G, h_shelf0, h_shelf1, ddraft)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf0 !< the previous thickness of the ice shelf [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_shelf1 !< the current thickness of the ice shelf [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout)    :: ddraft !< the change in shelf draft thickness
  real :: b0,b1
  integer :: i, j, isc, iec, jsc, jec
  real    :: rhoi_rhow, OD

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  ddraft = 0.0

  do j=jsc,jec
    do i=isc,iec

      b0 = 0.0 ; b1 = 0.0

      if (h_shelf0(i,j)>0.0) then
        OD = CS%bed_elev(i,j) - rhoi_rhow * h_shelf0(i,j)
        if (OD >= 0) then
          !floating
          b0 = rhoi_rhow * h_shelf0(i,j)
        else
          b0 = CS%bed_elev(i,j)
        endif
      endif

      if (h_shelf1(i,j)>0.0) then
        OD = CS%bed_elev(i,j) - rhoi_rhow * h_shelf1(i,j)
        if (OD >= 0) then
          !floating
          b1 = rhoi_rhow * h_shelf1(i,j)
        else
          b1 = CS%bed_elev(i,j)
        endif
      endif

      ddraft(i,j) = b1-b0
    enddo
  enddo
end subroutine change_in_draft

!> This subroutine calculates the gradients of bilinear basis elements that
!! that are centered at the vertices of the cell.  Values are calculated at
!! points of gaussian quadrature.
subroutine bilinear_shape_functions (X, Y, Phi, area)
  real, dimension(4),   intent(in)    :: X   !< The x-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(4),   intent(in)    :: Y   !< The y-positions of the vertices of the quadrilateral [L ~> m].
  real, dimension(8,4), intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                             !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real,                 intent(out)   :: area !< The quadrilateral cell area [L2 ~> m2].

! X and Y must be passed in the form
    !  3 - 4
    !  |   |
    !  1 - 2

! this subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell. values are calculated at
! points of gaussian quadrature. (in 1D: .5 * (1 +/- sqrt(1/3)) for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear
!
! This should be a one-off; once per nonlinear solve? once per lifetime?
! ... will all cells have the same shape and dimension?

  real, dimension(4) :: xquad, yquad ! [nondim]
  real :: a,b,c,d  ! Various lengths [L ~> m]
  real :: xexp, yexp ! [nondim]
  integer :: node, qpoint, xnode, ynode

  xquad(1:3:2) = .5 * (1-sqrt(1./3)) ; yquad(1:2) = .5 * (1-sqrt(1./3))
  xquad(2:4:2) = .5 * (1+sqrt(1./3)) ; yquad(3:4) = .5 * (1+sqrt(1./3))

  do qpoint=1,4

    a = ((-X(1)*(1-yquad(qpoint)))+(X(4)*yquad(qpoint))) + ((X(2)*(1-yquad(qpoint)))-(X(3)*yquad(qpoint))) !d(x)/d(x*)
    b = ((-Y(1)*(1-yquad(qpoint)))+(Y(4)*yquad(qpoint))) + ((Y(2)*(1-yquad(qpoint)))-(Y(3)*yquad(qpoint))) !d(y)/d(x*)
    c = ((-X(1)*(1-xquad(qpoint)))+(X(4)*xquad(qpoint))) + ((-X(2)*xquad(qpoint))+(X(3)*(1-xquad(qpoint))))!d(x)/d(y*)
    d = ((-Y(1)*(1-xquad(qpoint)))+(Y(4)*xquad(qpoint))) + ((-Y(2)*xquad(qpoint))+(Y(3)*(1-xquad(qpoint))))!d(y)/d(y*)

    do node=1,4

      xnode = 2-mod(node,2) ; ynode = ceiling(REAL(node)/2)

      if (ynode == 1) then
        yexp = 1-yquad(qpoint)
      else
        yexp = yquad(qpoint)
      endif

      if (1 == xnode) then
        xexp = 1-xquad(qpoint)
      else
        xexp = xquad(qpoint)
      endif

      Phi(2*node-1,qpoint) = ( d * (2 * xnode - 3) * yexp - b * (2 * ynode - 3) * xexp) / ((a*d)-(b*c))
      Phi(2*node,qpoint)   = (-c * (2 * xnode - 3) * yexp + a * (2 * ynode - 3) * xexp) / ((a*d)-(b*c))

    enddo
  enddo

  area = quad_area(X, Y)

end subroutine bilinear_shape_functions

!> This subroutine calculates the gradients of bilinear basis elements that are centered at the
!! vertices of the cell using a locally orthogoal MOM6 grid.  Values are calculated at
!! points of gaussian quadrature.
subroutine bilinear_shape_fn_grid(G, i, j, Phi, Jac)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  integer,               intent(in)    :: i   !< The i-index in the grid to work on.
  integer,               intent(in)    :: j   !< The j-index in the grid to work on.
  real, dimension(8,4),  intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                              !! quadrature points surrounding the cell vertices [L-1 ~> m-1].
  real, dimension(4), optional, intent(out) :: Jac !< Jacobian determinant |J_q| = a_q*d_q at each
                                                   !! Gaussian quadrature point [L2 ~> m2].

! This subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell.  The values are calculated at
! points of gaussian quadrature. (in 1D: .5 * (1 +/- sqrt(1/3)) for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1,j) gives d(Phi_i)/dx at quadrature point j
! Phi(2*i,j) gives d(Phi_i)/dy at quadrature point j
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear
!
! This should be a one-off; once per nonlinear solve? once per lifetime?

  real, dimension(4) :: xquad, yquad ! [nondim]
  ! Mirror lookups: xquad_m(qp) == 1 - xquad(qp), yquad_m(qp) == 1 - yquad(qp) mathematically,
  ! but each mirror entry is the stored value at the x- or y-mirrored quadrature point. This
  ! ensures rotation-paired QPs read bit-identical operand values (avoids the (1 - v) vs v_other
  ! 1-ulp asymmetry that breaks rotation invariance)
  real, dimension(4) :: xquad_m, yquad_m
  real :: a, d       ! Interpolated grid spacings [L ~> m]
  real :: xexp, yexp ! [nondim]
  integer :: node, qpoint, xnode, ynode

  xquad(1:3:2) = .5 * (1-sqrt(1./3)) ; yquad(1:2) = .5 * (1-sqrt(1./3))
  xquad(2:4:2) = .5 * (1+sqrt(1./3)) ; yquad(3:4) = .5 * (1+sqrt(1./3))

  ! x-mirror swaps qp 1<->2 and 3<->4; y-mirror swaps 1<->3 and 2<->4
  xquad_m(1) = xquad(2) ; xquad_m(2) = xquad(1) ; xquad_m(3) = xquad(4) ; xquad_m(4) = xquad(3)
  yquad_m(1) = yquad(3) ; yquad_m(2) = yquad(4) ; yquad_m(3) = yquad(1) ; yquad_m(4) = yquad(2)

  do qpoint=1,4
    ! Fallback uses the single available face length when at the global
    ! south/west edge (no neighbor to interpolate against) or when J-1 / I-1
    ! falls outside q-axis array bounds (nonsymmetric mode in deepest halo).
    if ((J-1 >= G%JsdB) .and. (j + G%jdg_offset > G%jsg)) then
      a = (G%dxCv(i,J-1) * yquad_m(qpoint)) + (G%dxCv(i,J) * yquad(qpoint)) ! d(x)/d(x*)
    else
      a = G%dxCv(i,J) !* yquad(qpoint) ! d(x)/d(x*)
    endif
    if ((I-1 >= G%IsdB) .and. (i + G%idg_offset > G%isg)) then
      d = (G%dyCu(I-1,j) * xquad_m(qpoint)) + (G%dyCu(I,j) * xquad(qpoint)) ! d(y)/d(y*)
    else
      d = G%dyCu(I,j) !* xquad(qpoint)
    endif

    do node=1,4
      xnode = 2-mod(node,2) ; ynode = ceiling(REAL(node)/2)

      if (ynode == 1) then
        yexp = yquad_m(qpoint)
      else
        yexp = yquad(qpoint)
      endif

      if (1 == xnode) then
        xexp = xquad_m(qpoint)
      else
        xexp = xquad(qpoint)
      endif

      Phi(2*node-1,qpoint) = ( (d * (2 * xnode - 3)) * yexp ) / (a*d)
      Phi(2*node,qpoint)   = ( (a * (2 * ynode - 3)) * xexp ) / (a*d)

    enddo
    if (present(Jac)) Jac(qpoint) = a * d
  enddo

end subroutine bilinear_shape_fn_grid

!> This subroutine calculates the gradients of bilinear basis elements that are centered at the
!! vertices of the cell using a locally orthogoal MOM6 grid.  Values are calculated at
!! a sinlge cell-centered quadrature point, which should match the grid cell h-point
subroutine bilinear_shape_fn_grid_1qp(G, i, j, Phi)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  integer,               intent(in)    :: i   !< The i-index in the grid to work on.
  integer,               intent(in)    :: j   !< The j-index in the grid to work on.
  real, dimension(8),    intent(inout) :: Phi !< The gradients of bilinear basis elements at Gaussian
                                              !! quadrature points surrounding the cell vertices [L-1 ~> m-1].

! This subroutine calculates the gradients of bilinear basis elements that
! that are centered at the vertices of the cell.  The values are calculated at
! a cell-cented point of gaussian quadrature. (in 1D: .5 for [0,1])
!     (ordered in same way as vertices)
!
! Phi(2*i-1) gives d(Phi_i)/dx at the quadrature point
! Phi(2*i) gives d(Phi_i)/dy at the quadrature point
! Phi_i is equal to 1 at vertex i, and 0 at vertex k /= i, and bilinear

  real :: a, d       ! Interpolated grid spacings [L ~> m]
  real :: xexp=0.5, yexp=0.5 ! [nondim]
  integer :: node, qpoint, xnode, ynode

    ! d(x)/d(x*). Fallback to the single available face length at the global
    ! southern edge or when J-1 falls outside the q-axis array bounds (JsdB),
    ! which can happen in the deepest south halo on a nonsymmetric grid.
    if ((J-1 >= G%JsdB) .and. (j + G%jdg_offset > G%jsg)) then
      a = 0.5 * (G%dxCv(i,J-1) + G%dxCv(i,J))
    else
      a = G%dxCv(i,J)
    endif

    ! d(y)/d(y*). Same logic as above on the western edge / IsdB bound.
    if ((I-1 >= G%IsdB) .and. (i + G%idg_offset > G%isg)) then
      d = 0.5 * (G%dyCu(I-1,j) + G%dyCu(I,j))
    else
      d = G%dyCu(I,j)
    endif

    do node=1,4
      xnode = 2-mod(node,2) ; ynode = ceiling(REAL(node)/2)
      Phi(2*node-1) = ( (d * (2 * xnode - 3)) * yexp ) / (a*d)
      Phi(2*node)   = ( (a * (2 * ynode - 3)) * xexp ) / (a*d)
    enddo
end subroutine bilinear_shape_fn_grid_1qp


subroutine bilinear_shape_functions_subgrid(Phisub, nsub)
  integer, intent(in)    :: nsub   !< The number of subgridscale quadrature locations in each direction
  real, dimension(2,2,nsub,nsub,2,2), &
           intent(inout) :: Phisub !< Quadrature structure weights at subgridscale
                                   !! locations for finite element calculations [nondim]

  ! this subroutine is a helper for interpolation of floatation condition
  ! for the purposes of evaluating the terms \int (u,v) \phi_i dx dy in a cell that is
  !     in partial floatation
  ! the array Phisub contains the values of \phi_i (where i is a node of the cell)
  !     at quad point j
  ! i think this general approach may not work for nonrectangular elements...
  !

  ! Phisub(q1,q2,i,j,k,l)
  !  q1: quad point x-index
  !  q2: quad point y-index
  !  i: subgrid index in x-direction
  !  j: subgrid index in y-direction
  !  k: basis function x-index
  !  l: basis function y-index

  ! e.g. k=1,l=1 => node 1
  !      q1=2,q2=1 => quad point 2

    !  3 - 4
    !  |   |
    !  1 - 2

  integer :: i, j, qx, qy
  real,dimension(2)    :: xquad ! [nondim]
  real                 :: fracx ! The fractional sub-cell area in reference space [nondim]
  ! Mirror-symmetric per-direction node weights: a_left == 1-x_global, a_right == x_global
  ! mathematically, but constructed so that a_right(qx,i) is computed by exactly the same
  ! operand sequence as a_left(3-qx, nsub+1-i). This guarantees bit-exact rotation symmetry.
  real, dimension(2,nsub) :: a_left, a_right ! [nondim]

  xquad(1) = .5 * (1-sqrt(1./3)) ; xquad(2) = .5 * (1+sqrt(1./3))
  fracx = 1.0/real(nsub)

  do i=1,nsub ; do qx=1,2
    a_left (qx,i) = (real(nsub-i) + xquad(3-qx)) * fracx
    a_right(qx,i) = (real(i-1)    + xquad(qx))   * fracx
  enddo ; enddo

  do j=1,nsub ; do i=1,nsub
    do qy=1,2 ; do qx=1,2
      Phisub(qx,qy,i,j,1,1) = a_left (qx,i) * a_left (qy,j)
      Phisub(qx,qy,i,j,1,2) = a_left (qx,i) * a_right(qy,j)
      Phisub(qx,qy,i,j,2,1) = a_right(qx,i) * a_left (qy,j)
      Phisub(qx,qy,i,j,2,2) = a_right(qx,i) * a_right(qy,j)
    enddo ; enddo
  enddo ; enddo

end subroutine bilinear_shape_functions_subgrid


subroutine update_velocity_masks(CS, G, hmask, umask, vmask, u_face_mask, v_face_mask)
  type(ice_shelf_dyn_CS),intent(in)    :: CS !< A pointer to the ice shelf dynamics control structure
  type(ocean_grid_type), intent(inout) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: umask !< A coded mask indicating the nature of the
                                             !! zonal flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: vmask !< A coded mask indicating the nature of the
                                             !! meridional flow at the corner point
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: u_face_mask !< A coded mask for velocities at the C-grid u-face
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(out)   :: v_face_mask !< A coded mask for velocities at the C-grid v-face
  ! sets masks for velocity solve
  ! ignores the fact that their might be ice-free cells - this only considers the computational boundary

  ! !!!IMPORTANT!!! relies on thickness mask - assumed that this is called after hmask has been updated & halo-updated

  integer :: i, j, k, iscq, iecq, jscq, jecq, isd, jsd, is, js, iegq, jegq
  integer :: giec, gjec, gisc, gjsc, isc, jsc, iec, jec

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  iscq = G%iscB ; iecq = G%iecB ; jscq = G%jscB ; jecq = G%jecB
  isd = G%isd ; jsd = G%jsd
  iegq = G%iegB ; jegq = G%jegB
  gisc = G%Domain%nihalo ; gjsc = G%Domain%njhalo
  giec = G%Domain%niglobal+gisc ; gjec = G%Domain%njglobal+gjsc

  umask(:,:) = 0 ; vmask(:,:) = 0
  u_face_mask(:,:) = 0 ; v_face_mask(:,:) = 0

  if (G%symmetric) then
    is = isd ; js = jsd
  else
    is = isd+1 ; js = jsd+1
  endif

  do j=js,G%jed ; do i=is,G%ied
    if (hmask(i,j) == 1 .or. hmask(i,j)==3) then
      umask(I-1:I,J-1:J)=1
      vmask(I-1:I,J-1:J)=1
    endif
  enddo ; enddo

  do j=js,G%jed
    do i=is,G%ied

      if ((hmask(i,j) == 1) .OR. (hmask(i,j) == 3)) then

        do k=0,1

          select case (int(CS%u_face_mask_bdry(I-1+k,j)))
            case (5)
              umask(I-1+k,J-1:J) = 3.
              u_face_mask(I-1+k,j) = 5.
            case (3)
              umask(I-1+k,J-1:J) = 3.
              vmask(I-1+k,J-1:J) = 3.
              u_face_mask(I-1+k,j) = 3.
            case (6)
              vmask(I-1+k,J-1:J) = 3.
              u_face_mask(I-1+k,j) = 6.
            case (2)
              u_face_mask(I-1+k,j) = 2.
            case (4)
              umask(I-1+k,J-1:J) = 0.
              u_face_mask(I-1+k,j) = 4.
            case (0)
              umask(I-1+k,J-1:J) = 0.
              u_face_mask(I-1+k,j) = 0.
            case (1)  ! stress free x-boundary
              umask(I-1+k,J-1:J) = 0.
            case default
              umask(I-1+k,J-1) = max(1. , umask(I-1+k,J-1))
              umask(I-1+k,J)   = max(1. , umask(I-1+k,J))
          end select
        enddo

        do k=0,1

          select case (int(CS%v_face_mask_bdry(i,J-1+k)))
            case (5)
              vmask(I-1:I,J-1+k) = 3.
              v_face_mask(i,J-1+k) = 5.
            case (3)
              vmask(I-1:I,J-1+k) = 3.
              umask(I-1:I,J-1+k) = 3.
              v_face_mask(i,J-1+k) = 3.
            case (6)
              umask(I-1:I,J-1+k) = 3.
              v_face_mask(i,J-1+k) = 6.
            case (2)
              v_face_mask(i,J-1+k) = 2.
            case (4)
              vmask(I-1:I,J-1+k) = 0.
              v_face_mask(i,J-1+k) = 4.
            case (0)
              vmask(I-1:I,J-1+k) = 0.
              v_face_mask(i,J-1+k) = 0.
            case (1) ! stress free y-boundary
              vmask(I-1:I,J-1+k) = 0.
            case default
              vmask(I-1,J-1+k) = max(1. , vmask(I-1,J-1+k))
              vmask(I,J-1+k)   = max(1. , vmask(I,J-1+k))
          end select
        enddo


        if (i < G%ied) then
          if ((hmask(i+1,j) == 0) .OR. (hmask(i+1,j) == 2)) then
            ! east boundary or adjacent to unfilled cell
            u_face_mask(I,j) = 2.
          endif
        endif

        if (i > G%isd) then
          if ((hmask(i-1,j) == 0) .OR. (hmask(i-1,j) == 2)) then
            !adjacent to unfilled cell
            u_face_mask(I-1,j) = 2.
          endif
        endif

        if (j > G%jsd) then
          if ((hmask(i,j-1) == 0) .OR. (hmask(i,j-1) == 2)) then
            !adjacent to unfilled cell
            v_face_mask(i,J-1) = 2.
          endif
        endif

        if (j < G%jed) then
          if ((hmask(i,j+1) == 0) .OR. (hmask(i,j+1) == 2)) then
            !adjacent to unfilled cell
            v_face_mask(i,j) = 2.
          endif
        endif


      endif

    enddo
  enddo

  ! note: if the grid is nonsymmetric, there is a part that will not be transferred with a halo update
  ! so this subroutine must update its own symmetric part of the halo

  call pass_vector(u_face_mask, v_face_mask, G%domain, TO_ALL, CGRID_NE)
  call pass_vector(umask, vmask, G%domain, TO_ALL, BGRID_NE)

end subroutine update_velocity_masks

!> Interpolate the ice shelf thickness from tracer point to nodal points,
!! subject to a mask.
subroutine interpolate_H_to_B(G, h_shelf, hmask, H_node, min_h_shelf)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: h_shelf !< The ice shelf thickness at tracer points [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: H_node !< The ice shelf thickness at nodal (corner)
                                             !! points [Z ~> m].
  real, intent(in) :: min_h_shelf !< The minimum ice thickness used during ice dynamics [Z ~> m].

  integer :: i, j, isc, iec, jsc, jec, num_h, k, l, ic, jc
  real    :: h_arr(2,2)

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  H_node(:,:) = 0.0

  ! H_node is node-centered; average over all cells that share that node
  ! if no (active) cells share the node then its value there is irrelevant

  do j=jsc-1,jec
    do i=isc-1,iec
      num_h = 0
      do l=1,2 ; jc=j-1+l ; do k=1,2 ; ic=i-1+k
        if (hmask(ic,jc) == 1.0 .or. hmask(ic,jc) == 3.0) then
          h_arr(k,l)=max(h_shelf(ic,jc),min_h_shelf)
          num_h = num_h + 1
        else
          h_arr(k,l)=0.0
        endif
        if (num_h > 0) then
          H_node(i,j) = ((h_arr(1,1)+h_arr(2,2))+(h_arr(1,2)+h_arr(2,1))) / num_h
        endif
      enddo ; enddo
    enddo
  enddo

  call pass_var(H_node, G%domain,position=CORNER)

end subroutine interpolate_H_to_B

!> DG(1)-aware variant of interpolate_H_to_B. Each B-grid corner is shared
!! by up to four T-cells; in each cell we evaluate the nodal Q1 thickness
!! at the corresponding reference corner and average over the cells that
!! contribute (using the same hmask + min_h_shelf floor rules as
!! interpolate_H_to_B).
subroutine interpolate_H_to_B_DG(G, h_shelf, h_nodal, hmask, H_node, min_h_shelf)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: h_shelf !< Area-weighted cell-mean ice thickness Hbar [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G),2,2), &
                         intent(in)    :: h_nodal !< Nodal Q1 thickness at 4 corners per cell [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: hmask !< Ice shelf mask.
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: H_node !< Ice shelf thickness at nodal (corner) points [Z ~> m].
  real, intent(in) :: min_h_shelf !< The minimum ice thickness used during ice dynamics [Z ~> m].

  integer :: i, j, isc, iec, jsc, jec, num_h, k, l, ic, jc
  real    :: h_arr(2,2)

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  H_node(:,:) = 0.0

  ! Nodal Q1 reconstruction at B-grid node (i,j). The node is shared by up
  ! to 4 T-cells; each cell contributes the value at its own corresponding
  ! corner (R5 / R25 of nodal plan):
  !   cell (i,  j  ) at its NE = h_nodal(i,  j,  2,2)
  !   cell (i+1,j  ) at its NW = h_nodal(i+1,j,  1,2)
  !   cell (i,  j+1) at its SE = h_nodal(i,  j+1,2,1)
  !   cell (i+1,j+1) at its SW = h_nodal(i+1,j+1,1,1)
  ! Loop indexing (k,l) maps to (ic=i-1+k, jc=j-1+l); see comment block in
  ! the original modal version for the reference-corner mapping.
  do j=jsc-1,jec
    do i=isc-1,iec
      num_h = 0
      h_arr(:,:) = 0.0
      do l=1,2 ; jc=j-1+l ; do k=1,2 ; ic=i-1+k
        if (hmask(ic,jc) == 1.0 .or. hmask(ic,jc) == 3.0) then
          ! For cell (ic,jc), pick the corner at reference (xi,eta) =
          !   (k=1,l=1) -> (+0.5,+0.5) -> NE -> h_nodal(ic,jc,2,2)
          !   (k=2,l=1) -> (-0.5,+0.5) -> NW -> h_nodal(ic,jc,1,2)
          !   (k=1,l=2) -> (+0.5,-0.5) -> SE -> h_nodal(ic,jc,2,1)
          !   (k=2,l=2) -> (-0.5,-0.5) -> SW -> h_nodal(ic,jc,1,1)
          if (k == 1 .and. l == 1) then
            h_arr(k,l) = max(h_nodal(ic,jc,2,2), min_h_shelf)
          elseif (k == 2 .and. l == 1) then
            h_arr(k,l) = max(h_nodal(ic,jc,1,2), min_h_shelf)
          elseif (k == 1 .and. l == 2) then
            h_arr(k,l) = max(h_nodal(ic,jc,2,1), min_h_shelf)
          else
            h_arr(k,l) = max(h_nodal(ic,jc,1,1), min_h_shelf)
          endif
          num_h = num_h + 1
        else
          h_arr(k,l) = 0.0
        endif
        if (num_h > 0) then
          H_node(i,j) = ((h_arr(1,1)+h_arr(2,2))+(h_arr(1,2)+h_arr(2,1))) / num_h
        endif
      enddo ; enddo
    enddo
  enddo

  call pass_var(H_node, G%domain, position=CORNER)

end subroutine interpolate_H_to_B_DG

!> Deallocates all memory associated with the ice shelf dynamics module
subroutine ice_shelf_dyn_end(CS)
  type(ice_shelf_dyn_CS), pointer   :: CS !< A pointer to the ice shelf dynamics control structure

  if (.not.associated(CS)) return

  deallocate(CS%u_shelf, CS%v_shelf)
  deallocate(CS%taudx_shelf, CS%taudy_shelf)
  deallocate(CS%sx_shelf, CS%sy_shelf)
  deallocate(CS%t_shelf, CS%tmask)
  deallocate(CS%u_bdry_val, CS%v_bdry_val)
  deallocate(CS%u_face_mask, CS%v_face_mask)
  deallocate(CS%u_flux_bdry_val, CS%v_flux_bdry_val)
  deallocate(CS%umask, CS%vmask)
  deallocate(CS%u_face_mask_bdry, CS%v_face_mask_bdry)
  deallocate(CS%h_bdry_val)
  if (associated(CS%calve_mask)) deallocate(CS%calve_mask)

  deallocate(CS%ice_visc, CS%AGlen_visc)
  deallocate(CS%newton_visc_factor, CS%newton_str_ux, CS%newton_str_vy, CS%newton_str_sh)
  deallocate(CS%newton_umid, CS%newton_vmid, CS%newton_drag_coef)
  deallocate(CS%C_basal_friction)
  deallocate(CS%coef_prefactor, CS%fB_elem)
  deallocate(CS%OD_rt, CS%OD_av)
  deallocate(CS%t_bdry_val, CS%bed_elev, CS%bed_node)
  if (associated(CS%h_nodal)) deallocate(CS%h_nodal)
  if (associated(CS%h_flot)) deallocate(CS%h_flot)
  if (associated(CS%Minv_xi)) deallocate(CS%Minv_xi)
  if (associated(CS%Minv_eta)) deallocate(CS%Minv_eta)
  if (associated(CS%cell_mean_w)) deallocate(CS%cell_mean_w)
  if (associated(CS%h_source_rate)) deallocate(CS%h_source_rate)
  if (associated(CS%h_source_rate_last)) deallocate(CS%h_source_rate_last)
  if (associated(CS%phi_x_FV)) deallocate(CS%phi_x_FV)
  if (associated(CS%phi_y_FV)) deallocate(CS%phi_y_FV)
  if (associated(CS%dg_art_visc_coef_u)) deallocate(CS%dg_art_visc_coef_u)
  if (associated(CS%dg_art_visc_coef_v)) deallocate(CS%dg_art_visc_coef_v)
  if (associated(CS%dg_art_visc_nu_u)) deallocate(CS%dg_art_visc_nu_u)
  if (associated(CS%dg_art_visc_nu_v)) deallocate(CS%dg_art_visc_nu_v)
  if (associated(CS%dg_art_visc_excess_frac_u)) deallocate(CS%dg_art_visc_excess_frac_u)
  if (associated(CS%dg_art_visc_excess_frac_v)) deallocate(CS%dg_art_visc_excess_frac_v)
  if (associated(CS%dg_art_visc_allow_u)) deallocate(CS%dg_art_visc_allow_u)
  if (associated(CS%dg_art_visc_allow_v)) deallocate(CS%dg_art_visc_allow_v)
  if (associated(CS%dg_art_visc_cell_scale)) deallocate(CS%dg_art_visc_cell_scale)
  if (associated(CS%dg_slow_idle_face_u)) deallocate(CS%dg_slow_idle_face_u)
  if (associated(CS%dg_slow_idle_face_v)) deallocate(CS%dg_slow_idle_face_v)
  if (allocated(CS%mu_lim_xi))    deallocate(CS%mu_lim_xi)
  if (allocated(CS%mu_lim_eta))   deallocate(CS%mu_lim_eta)
  if (allocated(CS%mu_lim_cross)) deallocate(CS%mu_lim_cross)
  if (associated(CS%dg_lim_phi_xi))     deallocate(CS%dg_lim_phi_xi)
  if (associated(CS%dg_lim_phi_eta))    deallocate(CS%dg_lim_phi_eta)
  if (associated(CS%dg_lim_phi_cross))  deallocate(CS%dg_lim_phi_cross)
  if (associated(CS%dg_lim_mass_drift)) deallocate(CS%dg_lim_mass_drift)
  if (associated(CS%dg_lim_phi))        deallocate(CS%dg_lim_phi)
  if (associated(CS%dg_lim_pk_factor))  deallocate(CS%dg_lim_pk_factor)
  deallocate(CS%ground_frac, CS%ground_frac_rt)
  if (associated(CS%Jac)) deallocate(CS%Jac)
  if (associated(CS%Phi)) deallocate(CS%Phi)
  if (associated(CS%Phisub)) deallocate(CS%Phisub)
  if (associated(CS%PhiC)) deallocate(CS%PhiC)

  deallocate(CS)

end subroutine ice_shelf_dyn_end


!> This subroutine updates the vertically averaged ice shelf temperature.
subroutine ice_shelf_temp(CS, ISS, G, US, time_step, melt_rate, Time)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< A pointer to the ice shelf control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: melt_rate !< basal melt rate [R Z T-1 ~> kg m-2 s-1]
  type(time_type),        intent(in)    :: Time !< The current model time

!    This subroutine takes the velocity (on the Bgrid) and timesteps
!      (HT)_t = - div (uHT) + (adot Tsurf -bdot Tbot) once and then calculates T=HT/H
!
!    The flux overflows are included here. That is because they will be used to advect 3D scalars
!    into partial cells

  real, dimension(SZDI_(G),SZDJ_(G))   :: th_after_uflux, th_after_vflux, TH ! Integrated temperatures [C Z ~> degC m]
  integer                           :: isd, ied, jsd, jed, i, j, isc, iec, jsc, jec
  real :: Tsurf ! Surface air temperature [C ~> degC].  This is hard coded but should be an input argument.
  real :: adot  ! A surface heat exchange coefficient [R Z T-1 ~> kg m-2 s-1].


  ! For now adot and Tsurf are defined here adot=surf acc 0.1m/yr, Tsurf=-20oC, vary them later
  adot = (0.1/(365.0*86400.0))*US%m_to_Z*US%T_to_s * CS%density_ice
  Tsurf = -20.0*US%degC_to_C

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  th_after_uflux(:,:) = 0.0
  th_after_vflux(:,:) = 0.0

  do j=jsd,jed ; do i=isd,ied
!    if (ISS%hmask(i,j) > 1) then
    if ((ISS%hmask(i,j) == 3) .or. (ISS%hmask(i,j) == -2)) then
      CS%t_shelf(i,j) = CS%t_bdry_val(i,j)
    endif
  enddo ; enddo

  do j=jsd,jed ; do i=isd,ied
    ! Convert the averge temperature to a depth integrated temperature.
    TH(i,j) = CS%t_shelf(i,j)*ISS%h_shelf(i,j)
  enddo ; enddo


  call ice_shelf_advect_temp_x(CS, G, time_step, ISS%hmask, TH, th_after_uflux)
  call ice_shelf_advect_temp_y(CS, G, time_step, ISS%hmask, th_after_uflux, th_after_vflux)

  do j=jsc,jec ; do i=isc,iec
    ! Convert the integrated temperature back to the average temperature.
!   if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
    if (ISS%h_shelf(i,j) > 0.0) then
      CS%t_shelf(i,j) = th_after_vflux(i,j) / ISS%h_shelf(i,j)
    else
      CS%t_shelf(i,j) = CS%T_shelf_missing
    endif
!   endif

    if ((ISS%hmask(i,j) == 1) .or. (ISS%hmask(i,j) == 2)) then
      if (ISS%h_shelf(i,j) > 0.0) then
        CS%t_shelf(i,j) = CS%t_shelf(i,j) + &
            time_step*(adot*Tsurf - melt_rate(i,j)*ISS%tfreeze(i,j))/(CS%density_ice*ISS%h_shelf(i,j))
      else
        ! the ice is about to melt away in this case set thickness, area, and mask to zero
        ! NOTE: not mass conservative, should maybe scale salt & heat flux for this cell
        CS%t_shelf(i,j) = CS%T_shelf_missing
        CS%tmask(i,j) = 0.0
      endif
    elseif (ISS%hmask(i,j) == 0) then
      CS%t_shelf(i,j) = CS%T_shelf_missing
    elseif ((ISS%hmask(i,j) == 3) .or. (ISS%hmask(i,j) == -2)) then
      CS%t_shelf(i,j) = CS%t_bdry_val(i,j)
    endif
  enddo ; enddo

  call pass_var(CS%t_shelf, G%domain, complete=.false.)
  call pass_var(CS%tmask, G%domain, complete=.true.)

  if (CS%debug) then
    call hchksum(CS%t_shelf, "temp after front", G%HI, haloshift=3, unscale=US%C_to_degC)
  endif

end subroutine ice_shelf_temp


subroutine ice_shelf_advect_temp_x(CS, G, time_step, hmask, h0, h_after_uflux)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h0 !< The initial ice shelf thicknesses times temperature [C Z ~> degC m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_uflux !< The ice shelf thicknesses times temperature after
                                              !! the zonal mass fluxes [C Z ~> degC m]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization

  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  integer :: i_off, j_off
  logical :: at_east_bdry, at_west_bdry
  real, dimension(-2:2) :: stencil ! A copy of the neighboring thicknesses times temperatures [C Z ~> degC m]
  real :: u_face     ! Zonal velocity at a face, positive if out [L T-1 ~> m s-1]
  real :: flux_diff  ! The difference in fluxes [C Z ~> degC m]
  real :: phi        ! A limiting ratio [nondim]

  is = G%isc-2 ; ie = G%iec+2 ; js = G%jsc ; je = G%jec ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do j=jsd+1,jed-1
    if (((j+j_off) <= G%domain%njglobal+G%domain%njhalo) .AND. &
        ((j+j_off) >= G%domain%njhalo+1)) then ! based on mehmet's code - only if btw north & south boundaries

      stencil(:) = 0.0 ! This is probably unnecessary, as the code is written
!     if (i+i_off == G%domain%nihalo+G%domain%nihalo)
      do i=is,ie

        if (((i+i_off) <= G%domain%niglobal+G%domain%nihalo) .AND. &
             ((i+i_off) >= G%domain%nihalo+1)) then

          if (i+i_off == G%domain%nihalo+1) then
            at_west_bdry=.true.
          else
            at_west_bdry=.false.
          endif

          if (i+i_off == G%domain%niglobal+G%domain%nihalo) then
            at_east_bdry=.true.
          else
            at_east_bdry=.false.
          endif

          if (hmask(i,j) == 1) then

            h_after_uflux(i,j) = h0(i,j)

            stencil(:) = h0(i-2:i+2,j)  ! fine as long has nx_halo >= 2

            flux_diff = 0

            ! 1ST DO LEFT FACE

            if (CS%u_face_mask(I-1,j) == 4.) then

              flux_diff = flux_diff + G%dyCu(I-1,j) * time_step * CS%u_flux_bdry_val(I-1,j) * &
                               CS%t_bdry_val(i-1,j) / G%areaT(i,j)
            else

              ! get u-velocity at center of left face
              u_face = 0.5 * (CS%u_shelf(I-1,J-1) + CS%u_shelf(I-1,J))

              if (u_face > 0) then !flux is into cell - we need info from h(i-2), h(i-1) if available

              ! i may not cover all the cases.. but i cover the realistic ones

                if (at_west_bdry .AND. (hmask(i-1,j) == 3)) then ! at western bdry but there is a
                              ! thickness bdry condition, and the stencil contains it
                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I-1,j) * time_step * stencil(-1) / G%areaT(i,j)

                elseif (hmask(i-1,j) * hmask(i-2,j) == 1) then  ! h(i-2) and h(i-1) are valid
                  phi = slope_limiter(stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff = flux_diff + ((ABS(u_face) * G%dyCu(I-1,j)* time_step / G%areaT(i,j)) * &
                           (stencil(-1) - (phi * (stencil(-1)-stencil(0))/2)))

                else                            ! h(i-1) is valid
                                    ! (o.w. flux would most likely be out of cell)
                                    !  but h(i-2) is not

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j) * stencil(-1)

                endif

              elseif (u_face < 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available
                if (hmask(i-1,j) * hmask(i+1,j) == 1) then         ! h(i-1) and h(i+1) are both valid
                  phi = slope_limiter(stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j)) * &
                             (stencil(0) - (phi * (stencil(0)-stencil(-1))/2)))

                else
                  flux_diff = flux_diff - ABS(u_face) * G%dyCu(I-1,j) * time_step / G%areaT(i,j) * stencil(0)
                endif
              endif
            endif

            ! NEXT DO RIGHT FACE

            ! get u-velocity at center of eastern face

            if (CS%u_face_mask(I,j) == 4.) then

              flux_diff = flux_diff + G%dyCu(I,j) * time_step * CS%u_flux_bdry_val(I,j) *&
                               CS%t_bdry_val(i+1,j) / G%areaT(i,j)
            else

              u_face = 0.5 * (CS%u_shelf(I,J-1) + CS%u_shelf(I,J))

              if (u_face < 0) then !flux is into cell - we need info from h(i+2), h(i+1) if available

                if (at_east_bdry .AND. (hmask(i+1,j) == 3)) then ! at eastern bdry but there is a
                                            ! thickness bdry condition, and the stencil contains it

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I,j) * time_step * stencil(1) / G%areaT(i,j)

                elseif (hmask(i+1,j) * hmask(i+2,j) == 1) then  ! h(i+2) and h(i+1) are valid

                  phi = slope_limiter(stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff = flux_diff + ((ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j)) * &
                      (stencil(1) - (phi * (stencil(1)-stencil(0))/2)))

                else                            ! h(i+1) is valid
                                            ! (o.w. flux would most likely be out of cell)
                                            !  but h(i+2) is not

                  flux_diff = flux_diff + ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j) * stencil(1)

                endif

              elseif (u_face > 0) then !flux is out of cell - we need info from h(i-1), h(i+1) if available

                if (hmask(i-1,j) * hmask(i+1,j) == 1) then         ! h(i-1) and h(i+1) are both valid

                  phi = slope_limiter(stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j)) * &
                      (stencil(0) - (phi * (stencil(0)-stencil(1))/2)))

                else  ! h(i+1) is valid (o.w. flux would most likely be out of cell) but h(i+2) is not

                  flux_diff = flux_diff - ABS(u_face) * G%dyCu(I,j) * time_step / G%areaT(i,j) * stencil(0)

                endif

              endif

              h_after_uflux(i,j) = h_after_uflux(i,j) + flux_diff

            endif

          endif

        endif

      enddo ! i loop

    endif

  enddo ! j loop

end subroutine ice_shelf_advect_temp_x

subroutine ice_shelf_advect_temp_y(CS, G, time_step, hmask, h_after_uflux, h_after_vflux)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< A pointer to the ice shelf control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_after_uflux !< The ice shelf thicknesses times temperature after
                                              !! the zonal mass fluxes [C Z ~> degC m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_after_vflux !< The ice shelf thicknesses times temperature after
                                              !! the meridional mass fluxes [C Z ~> degC m]

  ! use will be made of ISS%hmask here - its value at the boundary will be zero, just like uncovered cells
  ! if there is an input bdry condition, the thickness there will be set in initialization

  integer :: i, j, is, ie, js, je, isd, ied, jsd, jed
  integer :: i_off, j_off
  logical :: at_north_bdry, at_south_bdry
  real, dimension(-2:2) :: stencil ! A copy of the neighboring thicknesses times temperatures [C Z ~> degC m]
  real :: v_face     ! Pseudo-meridional velocity at a cell face, positive if out [L T-1 ~> m s-1]
  real :: flux_diff  ! The difference in fluxes [C Z ~> degC m]
  real :: phi

  is = G%isc ; ie = G%iec ; js = G%jsc-1 ; je = G%jec+1 ; isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset

  do i=isd+2,ied-2
    if (((i+i_off) <= G%domain%niglobal+G%domain%nihalo) .AND. &
       ((i+i_off) >= G%domain%nihalo+1)) then  ! based on mehmet's code - only if btw east & west boundaries

      stencil(:) = 0.0 ! This is probably unnecessary, as the code is written

      do j=js,je

        if (((j+j_off) <= G%domain%njglobal+G%domain%njhalo) .AND. &
             ((j+j_off) >= G%domain%njhalo+1)) then

          if (j+j_off == G%domain%njhalo+1) then
            at_south_bdry=.true.
          else
            at_south_bdry=.false.
          endif
          if (j+j_off == G%domain%njglobal+G%domain%njhalo) then
            at_north_bdry=.true.
          else
            at_north_bdry=.false.
          endif

          if (hmask(i,j) == 1) then
            h_after_vflux(i,j) = h_after_uflux(i,j)

            stencil(:) = h_after_uflux(i,j-2:j+2)  ! fine as long has ny_halo >= 2
            flux_diff = 0

            ! 1ST DO south FACE

            if (CS%v_face_mask(i,J-1) == 4.) then

              flux_diff = flux_diff + G%dxCv(i,J-1) * time_step * CS%v_flux_bdry_val(i,J-1) * &
                                 CS%t_bdry_val(i,j-1)/ G%areaT(i,j)
            else

              ! get u-velocity at center of west face
              v_face = 0.5 * (CS%v_shelf(I-1,J-1) + CS%v_shelf(I,J-1))

              if (v_face > 0) then !flux is into cell - we need info from h(j-2), h(j-1) if available

                ! i may not cover all the cases.. but i cover the realistic ones

                if (at_south_bdry .AND. (hmask(i,j-1) == 3)) then ! at western bdry but there is a
                                            ! thickness bdry condition, and the stencil contains it
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J-1) * time_step * stencil(-1) / G%areaT(i,j)

                elseif (hmask(i,j-1) * hmask(i,j-2) == 1) then  ! h(j-2) and h(j-1) are valid

                  phi = slope_limiter(stencil(-1)-stencil(-2), stencil(0)-stencil(-1))
                  flux_diff = flux_diff + ((ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j)) * &
                      (stencil(-1) - (phi * (stencil(-1)-stencil(0))/2)))

                else     ! h(j-1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j-2) is not
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j) * stencil(-1)
                endif

              elseif (v_face < 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) == 1) then  ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter(stencil(0)-stencil(1), stencil(-1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j)) * &
                      (stencil(0) - (phi * (stencil(0)-stencil(-1))/2)))
                else
                  flux_diff = flux_diff - ABS(v_face) * G%dxCv(i,J-1) * time_step / G%areaT(i,j) * stencil(0)
                endif

              endif

            endif

            ! NEXT DO north FACE

            if (CS%v_face_mask(i,J) == 4.) then
              flux_diff = flux_diff + G%dxCv(i,J) * time_step * CS%v_flux_bdry_val(i,J) *&
                               CS%t_bdry_val(i,j+1)/ G%areaT(i,j)
            else

            ! get u-velocity at center of east face
              v_face = 0.5 * (CS%v_shelf(I-1,J) + CS%v_shelf(I,J))

              if (v_face < 0) then !flux is into cell - we need info from h(j+2), h(j+1) if available

                if (at_north_bdry .AND. (hmask(i,j+1) == 3)) then ! at eastern bdry but there is a
                                            ! thickness bdry condition, and the stencil contains it
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J) * time_step * stencil(1) / G%areaT(i,j)
                elseif (hmask(i,j+1) * hmask(i,j+2) == 1) then  ! h(j+2) and h(j+1) are valid
                  phi = slope_limiter (stencil(1)-stencil(2), stencil(0)-stencil(1))
                  flux_diff = flux_diff + ((ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j)) * &
                      (stencil(1) - (phi * (stencil(1)-stencil(0))/2)))
                else     ! h(j+1) is valid
                         ! (o.w. flux would most likely be out of cell)
                         !  but h(j+2) is not
                  flux_diff = flux_diff + ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j) * stencil(1)
                endif

              elseif (v_face > 0) then !flux is out of cell - we need info from h(j-1), h(j+1) if available

                if (hmask(i,j-1) * hmask(i,j+1) == 1) then         ! h(j-1) and h(j+1) are both valid
                  phi = slope_limiter (stencil(0)-stencil(-1), stencil(1)-stencil(0))
                  flux_diff = flux_diff - ((ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j)) * &
                      (stencil(0) - (phi * (stencil(0)-stencil(1))/2)))
                else   ! h(j+1) is valid
                       ! (o.w. flux would most likely be out of cell)
                       !  but h(j+2) is not
                  flux_diff = flux_diff - ABS(v_face) * G%dxCv(i,J) * time_step / G%areaT(i,j) * stencil(0)
                endif

              endif

            endif

            h_after_vflux(i,j) = h_after_vflux(i,j) + flux_diff
          endif
        endif
      enddo ! j loop
    endif
  enddo ! i loop

end subroutine ice_shelf_advect_temp_y



!> Compute driving stress at B-grid nodes using a Pure DG(1) formulation.
!! Allows sub-element driving stress around grounding line. Evaluates the FEM
!! weak-form integral using integration by parts. Interior faces use a simple
!! central numerical flux P* = 0.5*(P_loc + P_ngh); any single-valued P* gives
!! the same SSA node-assembly total by the IBP-reverse identity, so penalty
!! terms (Rusanov, IIPG) were dropped. Dirichlet thickness BC sides use their
!! own P as authoritative.
subroutine calc_shelf_driving_stress_DG(CS, ISS, G, US, taudx, taudy, OD)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: taudx  !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: taudy  !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: OD  !< Ocean floor depth at tracer points [Z ~> m].

  ! Local variables
  real :: rho        ! Ice density [R ~> kg m-3]
  real :: rhow       ! Reference ocean density [R ~> kg m-3]
  real :: rhoi_rhow  ! Ice/ocean density ratio [nondim]
  real :: grav       ! Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]

  ! Face/Boundary variables
  real :: h_loc_A, h_loc_B ! Local thickness at face endpoints A and B [Z ~> m]
  real :: b_loc_A, b_loc_B ! Local bed elevation at face endpoints A and B [Z ~> m]
  real :: h_ngh_A, h_ngh_B ! Neighbor thickness at face endpoints A and B [Z ~> m]
  real :: b_ngh_A, b_ngh_B ! Neighbor bed elevation at face endpoints A and B [Z ~> m]
  real :: h_loc, b_loc     ! Quadrature-point local thickness and bed [Z ~> m]
  real :: h_ngh, b_ngh     ! Quadrature-point neighbor thickness and bed [Z ~> m]
  real :: P_loc, P_ngh     ! Local and Neighbor pressures [R Z L2 T-2 ~> kg s-2]
  real :: P_star           ! Numerical flux at face [R Z L2 T-2 ~> kg s-2]
  real :: d_ocean          ! Draft of the ice for ocean pressure calculation [Z ~> m]
  real :: t_face           ! Face quadrature parameter in [0,1] [nondim]
  real :: phi_A, phi_B     ! Linear nodal basis values at the face quadrature point [nondim]
  integer :: gp_face       ! Face Gauss-point loop index [nondim]
  logical :: is_ext_bdry   ! True if the face is an external (ocean) boundary
  real :: phi_val           ! Bilinear nodal basis value at a qp [nondim]
  ! Main-grid 2x2 Gauss quadrature variables
  real :: xi_gp, eta_gp     ! DG reference coordinates in [-0.5,0.5] at a main qp [nondim]
  real :: h_gp              ! Ice thickness at a qp [Z ~> m]
  real :: dhdx_gp, dhdy_gp  ! Thickness gradients at a qp, in physical coordinates [Z L-1 ~> nondim]
  real :: bed_gp            ! Bed elevation at a qp [Z ~> m]
  real :: dbdx_gp, dbdy_gp  ! Bed gradients at a qp, physical coords [Z L-1 ~> nondim]
  real :: dbdx_ref, dbdy_ref ! Bed gradients in reference coords (pre-Jacobian) [Z ~> m]
  real :: bottom_force_x, bottom_force_y ! Bed/water bottom drag forces [R Z L2 T-2 ~> kg s-2]
  real :: dphi_dx_ref, dphi_dy_ref ! Basis function derivatives in reference coordinates [nondim]
  real :: dphi_dx, dphi_dy  ! Basis function derivatives in physical coordinates [L-1 ~> m-1]
  real :: p_term_vol        ! Integrated-by-parts volume pressure term [R Z L2 T-2 ~> kg s-2]
  real :: a_qp, d_qp        ! Per-qp interpolated cell-edge spacings [L ~> m]
  real :: weight            ! Per-qp quadrature weight including Jacobian [L2 ~> m2]
  real :: bed_corners(2,2)  ! Bed elevation at the 4 B-grid corners of an element [Z ~> m]
  real :: dxCv_S, dxCv_N    ! Cell-edge spacings on south and north faces [L ~> m]
  real :: dyCu_W, dyCu_E    ! Cell-edge spacings on west and east faces [L ~> m]

  ! Gauss quadrature: 2 nodes of a 2-point Gauss-Legendre rule on [0,1]
  real, dimension(2) :: xquad  ! Quadrature point locations on [0,1] [nondim]
  integer :: i, j, iq, jq, isc, iec, jsc, jec, isd, ied, jsd, jed
  integer :: i_off, j_off, gisc, giec, gjsc, gjec

  real, dimension(2,2,2,2) :: qp_dx, qp_dy ! Per-QP x and y stress contributions to the
                                           ! 4 cell-corners, indexed (qx,qy,m,n)
                                           ! [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2) :: vol_dx, vol_dy   ! Per-corner volume-integral driving-stress
                                           ! total for this cell, indexed (m,n)
                                           ! [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2) :: cell_dx_node, cell_dy_node ! Per-corner total (volume +
                                           ! Neumann face) for this cell, indexed (m,n)
                                           ! [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2) :: slope_x_gp, slope_y_gp ! Per-QP surface slopes pre-multiplied by
                                                 ! Jacobian a*d for physical-area weighting
                                                 ! [Z L ~> m]
  real, dimension(2,2) :: weight_gp ! Per-QP Jacobian a*d used for area-weighting slopes [L2 ~> m2]
  real :: slope_w_sum ! Sum of per-QP Jacobian weights for slope average [L2 ~> m2]
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: taudx_b, taudy_b
                                           !< Per-node 4-slot driving-stress
                                           !! accumulator, slot k indexed by which
                                           !! of the 4 surrounding cells contributes
                                           !! (1=SW, 2=SE, 3=NW, 4=NE)
                                           !! [R L3 Z T-2 ~> kg m s-2]

  ! Per-face Neumann-BC contributions to the face's two corner nodes (A,B).
  real :: face_dx_W_A, face_dx_W_B  ! West face x-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]
  real :: face_dx_E_A, face_dx_E_B  ! East face x-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]
  real :: face_dy_S_A, face_dy_S_B  ! South face y-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]
  real :: face_dy_N_A, face_dy_N_B  ! North face y-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]
  integer :: m, n  ! Cell-corner (m,n) loop indices in [1..2] [nondim]
  logical :: calc_slope_diag ! True if slope diagnostics will be calculated
  logical :: loc_is_bc       ! True if local cell has hmask==3 (Dirichlet thickness BC)
  logical :: ngh_is_bc       ! True if neighbor cell has hmask==3 (Dirichlet thickness BC)
  logical :: is_wall         ! True if the face is a velocity-Dirichlet wall (non-reentrant
                             ! global boundary, or explicit face_mask_bdry zero-normal-velocity
                             ! code), not flagged as a Neumann BC
  logical :: is_grounded     ! True if the local QP / face point is grounded (vs floating)

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset
  gisc = 1 ; gjsc = 1
  giec = G%domain%niglobal ; gjec = G%domain%njglobal

  rho = CS%density_ice
  rhow = CS%density_ocean_avg
  grav = CS%g_Earth
  rhoi_rhow = rho / rhow

  ! Gauss quadrature points on [0,1] (2-point rule)
  xquad(1) = 0.5 * (1.0 - sqrt(1.0/3.0))
  xquad(2) = 0.5 * (1.0 + sqrt(1.0/3.0))

  taudx(:,:) = 0.0 ; taudy(:,:) = 0.0
  taudx_b(:,:,:) = 0.0 ; taudy_b(:,:,:) = 0.0

  if (CS%id_sx_shelf > 0 .or.  CS%id_sy_shelf > 0 .or. CS%id_surf_slope_mag_shelf > 0) then
    calc_slope_diag=.true.
  else
    calc_slope_diag=.false.
  endif

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then

      ! Gather bed_node corners for this element
      ! Node ordering:  3(I-1,J) - 4(I,J)
      !                 |           |
      !                 1(I-1,J-1) - 2(I,J-1)
      bed_corners(1,1) = CS%bed_node(I-1,J-1)
      bed_corners(2,1) = CS%bed_node(I,J-1)
      bed_corners(1,2) = CS%bed_node(I-1,J)
      bed_corners(2,2) = CS%bed_node(I,J)

      ! Cell-edge metrics for Jacobian interpolation
      dxCv_S = G%dxCv(i,J-1) ; dxCv_N = G%dxCv(i,J)
      dyCu_W = G%dyCu(I-1,j) ; dyCu_E = G%dyCu(I,j)

      ! GL band: defer to sub-grid routine (mirrors CG_action_subgrid_basal).
      ! Non-GL cells: inline 2x2 main-grid Gauss quadrature, mirroring how
      ! CG_action handles the main-grid 4-qp loop. Both paths produce vol_d*(m,n).
      if (CS%GL_regularize .and. CS%ground_frac(i,j) > 0.0 .and. CS%ground_frac(i,j) < 1.0) then
        call calc_shelf_driving_stress_DG_subgrid(CS, CS%Phisub, &
            ISS%h_shelf(i,j), CS%h_nodal(i,j,:,:), bed_corners, &
            dxCv_S, dxCv_N, dyCu_W, dyCu_E, &
            rho, rhow, rhoi_rhow, grav, vol_dx, vol_dy, CS%sx_shelf(i,j), CS%sy_shelf(i,j), calc_slope_diag)
      else
        qp_dx(:,:,:,:) = 0.0 ; qp_dy(:,:,:,:) = 0.0
        do jq=1,2 ; do iq=1,2
          xi_gp  = xquad(iq) - 0.5
          eta_gp = xquad(jq) - 0.5

          a_qp = (dxCv_S * xquad(3-jq)) + (dxCv_N * xquad(jq))
          d_qp = (dyCu_W * xquad(3-iq)) + (dyCu_E * xquad(iq))
          weight = 0.25 * (a_qp * d_qp)

          ! Nodal Q1 evaluation of h and grad(h) at the QP (R9/R1 of nodal plan).
          ! Rotation-paired sum (SW+NE) + (SE+NW). xquad(iq), xquad(jq) in [0,1].
          h_gp = ((CS%h_nodal(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                  (CS%h_nodal(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
                 ((CS%h_nodal(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                  (CS%h_nodal(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
          h_gp = max(h_gp, CS%min_h_shelf)
          ! dN(a,b)/dxi at QP: dN(1,*)/dxi = -(eta or 1-eta); dN(2,*)/dxi = +.
          ! dN(a,b)/deta at QP: dN(*,1)/deta = -(xi or 1-xi); dN(*,2)/deta = +.
          dhdx_gp = ( ((-xquad(3-jq))*CS%h_nodal(i,j,1,1) + ( xquad(jq))   *CS%h_nodal(i,j,2,2)) + &
                      (( xquad(3-jq))*CS%h_nodal(i,j,2,1) + (-xquad(jq))   *CS%h_nodal(i,j,1,2)) ) / a_qp
          dhdy_gp = ( ((-xquad(3-iq))*CS%h_nodal(i,j,1,1) + ( xquad(iq))   *CS%h_nodal(i,j,2,2)) + &
                      ((-xquad(iq))  *CS%h_nodal(i,j,2,1) + ( xquad(3-iq)) *CS%h_nodal(i,j,1,2)) ) / d_qp

          bed_gp = ((bed_corners(1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                    (bed_corners(2,2) * (xquad(iq)   * xquad(jq))))  + &
                   ((bed_corners(2,1) * (xquad(iq)   * xquad(3-jq))) + &
                    (bed_corners(1,2) * (xquad(3-iq) * xquad(jq))))

          ! Reference-coord bed gradients
          dbdx_ref = ((bed_corners(1,1) * (-xquad(3-jq))) + &
                      (bed_corners(2,2) * ( xquad(jq))))  + &
                     ((bed_corners(2,1) * ( xquad(3-jq))) + &
                      (bed_corners(1,2) * (-xquad(jq))))
          dbdy_ref = ((bed_corners(1,1) * (-xquad(3-iq))) + &
                      (bed_corners(2,2) * ( xquad(iq))))  + &
                     ((bed_corners(2,1) * (-xquad(iq)))   + &
                      (bed_corners(1,2) * ( xquad(3-iq))))
          dbdx_gp = dbdx_ref / a_qp
          dbdy_gp = dbdy_ref / d_qp

          if (CS%GL_couple) then
            is_grounded = (CS%ground_frac(i,j) >= 1.0)
          else
            is_grounded = (rhoi_rhow * h_gp - bed_gp > 0.0)
          endif

          ! Flotation-branched IBP decomposition of -rho*g*h*grad(s):
          !   grounded: P = 0.5*rho*g*h^2, bottom_force = rho*g*h*grad(bed).
          !   floating: P = 0.5*(1 - rhoi_rhow)*rho*g*h^2, bottom_force = 0.
          ! This makes the IBP face term (volume IBP -> face) carry the
          ! correct (1 - rhoi_rhow) scaling at floating cells so the SSA node
          ! assembly matches the strong form's grad(s) decomposition. The
          ! previous full-P decomposition over-coupled h-jumps at floating
          ! interior faces by a factor of 1/(1 - rhoi_rhow).
          if (is_grounded) then
            p_term_vol = 0.5 * rho * grav * h_gp**2
            bottom_force_x = rho * grav * h_gp * dbdx_gp
            bottom_force_y = rho * grav * h_gp * dbdy_gp
          else
            p_term_vol = 0.5 * (1.0 - rhoi_rhow) * rho * grav * h_gp**2
            bottom_force_x = 0.0
            bottom_force_y = 0.0
          endif

          ! For slope diagnostics. Slope values are pre-multiplied by the per-QP
          ! Jacobian a_qp*d_qp so the QP sum below is the physical-area integral;
          ! we divide by the sum of those Jacobians to get the area-weighted mean.
          if (calc_slope_diag) then
            weight_gp(iq,jq) = a_qp * d_qp
            if (CS%GL_couple) then
              if (CS%ground_frac(i,j)<1) then
                slope_x_gp(iq,jq) = (1.0 - rhoi_rhow) * dhdx_gp * weight_gp(iq,jq)
                slope_y_gp(iq,jq) = (1.0 - rhoi_rhow) * dhdy_gp * weight_gp(iq,jq)
              else
                slope_x_gp(iq,jq) = (dhdx_gp - dbdx_gp) * weight_gp(iq,jq)
                slope_y_gp(iq,jq) = (dhdy_gp - dbdy_gp) * weight_gp(iq,jq)
              endif
            else
              if (rhoi_rhow * h_gp - bed_gp <= 0.0) then
                slope_x_gp(iq,jq) = (1.0 - rhoi_rhow) * dhdx_gp * weight_gp(iq,jq)
                slope_y_gp(iq,jq) = (1.0 - rhoi_rhow) * dhdy_gp * weight_gp(iq,jq)
              else
                slope_x_gp(iq,jq) = (dhdx_gp - dbdx_gp) * weight_gp(iq,jq)
                slope_y_gp(iq,jq) = (dhdy_gp - dbdy_gp) * weight_gp(iq,jq)
              endif
            endif
          endif

          ! Weak-form Volume Integration by Parts; p_term_vol set above
          ! by the flotation branch.
          do n=1,2 ; do m=1,2
            phi_val = (merge(xquad(iq), xquad(3-iq), m == 2)) * &
                      (merge(xquad(jq), xquad(3-jq), n == 2))
            dphi_dx_ref = (merge(1.0, -1.0, m == 2)) * &
                          (merge(xquad(jq), xquad(3-jq), n == 2))
            dphi_dy_ref = (merge(xquad(iq), xquad(3-iq), m == 2)) * &
                          (merge(1.0, -1.0, n == 2))
            dphi_dx = dphi_dx_ref / a_qp
            dphi_dy = dphi_dy_ref / d_qp

            ! Geometric correction for the IBP pressure term on non-rectangular (e.g. lat/lon)
            ! elements. On a bilinear element, a(eta) = dxCv_S*(1-eta) + dxCv_N*eta varies with
            ! eta, so the reference-space IBP of P*d(phi)/deta requires an extra term
            ! P*phi*da/deta to satisfy the discrete divergence theorem. Without it, constant P
            ! gives a spurious metric-dependent driving stress on lat/lon grids. The 0.25 factor
            ! equals weight/(a_qp*d_qp), converting the reference-space integral to the nodal sum.
            qp_dx(iq,jq,m,n) = (weight * ((dphi_dx * p_term_vol) + (phi_val * bottom_force_x))) &
                              + (0.25 * phi_val * p_term_vol * (dyCu_E - dyCu_W))
            qp_dy(iq,jq,m,n) = (weight * ((dphi_dy * p_term_vol) + (phi_val * bottom_force_y))) &
                              + (0.25 * phi_val * p_term_vol * (dxCv_N - dxCv_S))
          enddo ; enddo

        enddo ; enddo

        do n=1,2 ; do m=1,2
          vol_dx(m,n) = (qp_dx(1,1,m,n) + qp_dx(2,2,m,n)) + (qp_dx(1,2,m,n) + qp_dx(2,1,m,n))
          vol_dy(m,n) = (qp_dy(1,1,m,n) + qp_dy(2,2,m,n)) + (qp_dy(1,2,m,n) + qp_dy(2,1,m,n))
        enddo ; enddo
      endif

      if (calc_slope_diag) then
        if (.not. (CS%GL_regularize .and. CS%ground_frac(i,j) > 0.0 .and. CS%ground_frac(i,j) < 1.0)) then
          ! Physical-area-weighted mean: slope_*_gp already premultiplied by a*d
          slope_w_sum = (weight_gp(1,1)+weight_gp(2,2)) + (weight_gp(1,2)+weight_gp(2,1))
          CS%sx_shelf(i,j) = ((slope_x_gp(1,1)+slope_x_gp(2,2)) + (slope_x_gp(1,2)+slope_x_gp(2,1))) / slope_w_sum
          CS%sy_shelf(i,j) = ((slope_y_gp(1,1)+slope_y_gp(2,2)) + (slope_y_gp(1,2)+slope_y_gp(2,1))) / slope_w_sum
        endif
      endif

      face_dx_W_A = 0.0 ; face_dx_W_B = 0.0
      face_dx_E_A = 0.0 ; face_dx_E_B = 0.0
      face_dy_S_A = 0.0 ; face_dy_S_B = 0.0
      face_dy_N_A = 0.0 ; face_dy_N_B = 0.0

      ! ======================================================================
      ! West Face (I-1) of cell (i,j)
      ! Local cell is right (i), Neighbor cell is left (i-1)
      ! ======================================================================
      h_loc_A = max(CS%h_nodal(i,j,1,1), CS%min_h_shelf)
      h_loc_B = max(CS%h_nodal(i,j,1,2), CS%min_h_shelf)
      b_loc_A = bed_corners(1,1) ; b_loc_B = bed_corners(1,2)

      h_ngh_A = max(CS%h_nodal(i-1,j,2,1), CS%min_h_shelf)
      h_ngh_B = max(CS%h_nodal(i-1,j,2,2), CS%min_h_shelf)
      b_ngh_A = bed_corners(1,1) ; b_ngh_B = bed_corners(1,2)

      ! Dirichlet thickness BC: override face thicknesses with h_bdry_val on any hmask==3 side.
      loc_is_bc = (ISS%hmask(i,j) == 3)
      ngh_is_bc = (ISS%hmask(i-1,j) == 3)
      if (loc_is_bc) then
        h_loc_A = max(CS%h_bdry_val(i,j), CS%min_h_shelf) ; h_loc_B = h_loc_A
      endif
      if (ngh_is_bc) then
        h_ngh_A = max(CS%h_bdry_val(i-1,j), CS%min_h_shelf) ; h_ngh_B = h_ngh_A
      endif

      is_ext_bdry = ((CS%u_face_mask_bdry(I-1,j) == 2) .or. &
                    ((ISS%hmask(i-1,j) == 0 .or. ISS%hmask(i-1,j) == 2) .and. &
                     (CS%reentrant_x .or. (i+i_off /= gisc))))
      is_wall = (.not. is_ext_bdry) .and. ( &
                  ((.not. CS%reentrant_x) .and. (i+i_off == gisc)) .or. &
                  (int(CS%u_face_mask_bdry(I-1,j)) == 3) .or. &
                  (int(CS%u_face_mask_bdry(I-1,j)) == 5) )

      do gp_face=1,2
        t_face = xquad(gp_face)
        h_loc = (1.0 - t_face)*h_loc_A + t_face*h_loc_B
        b_loc = (1.0 - t_face)*b_loc_A + t_face*b_loc_B

        if (rhoi_rhow * h_loc - b_loc > 0.0) then
          P_loc = 0.5 * grav * rho * h_loc**2
        else
          P_loc = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_loc**2
        endif

        if (is_ext_bdry) then
          ! Ice-front Neumann (natural BC: (h*sigma_dev)*n = P_ice - P_ocean).
          ! With the flotation-branched volume IBP p_term_vol, the per-cell
          ! face contribution sums to (P_loc - P_star)*n*phi*dS, so:
          !   grounded: P_loc = 0.5*rho*g*h^2 and we set P_star = P_ocean
          !             to recover net (P_ice - P_ocean).
          !   floating: P_loc = (1 - rhoi_rhow)*0.5*rho*g*h^2 already
          !             equals (P_ice - P_ocean) for d_ocean = rhoi_rhow*h,
          !             so the IBP face P_star must vanish.
          if (rhoi_rhow * h_loc - b_loc > 0.0) then
            d_ocean = max(0.0, min(b_loc, rhoi_rhow * h_loc))
            P_star = 0.5 * grav * rhow * d_ocean**2
          else
            P_star = 0.0
          endif
        else if (is_wall) then
          ! Velocity-Dirichlet global wall: mirror p_term_vol so the IBP
          ! face term exactly cancels the local volume IBP face term. The
          ! flotation-branched P_loc equals p_term_vol per QP, so use it
          ! directly.
          P_star = P_loc
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Interior face: central P_star. Any single-valued P_star gives the
          ! same SSA node-assembly total; penalty/upwind terms cancel by
          ! the IBP-reverse identity. Dirichlet sides override with their own P.
          if (loc_is_bc .and. .not. ngh_is_bc) then
            P_star = P_loc
          elseif (ngh_is_bc .and. .not. loc_is_bc) then
            P_star = P_ngh
          else
            P_star = 0.5 * (P_loc + P_ngh)
          endif
        endif

        phi_A = 1.0 - t_face ; phi_B = t_face
        ! Face normal n_x = -1 -> Boundary Integral = + P_star
        face_dx_W_A = face_dx_W_A + 0.5 * G%dyCu(I-1,j) * phi_A * P_star
        face_dx_W_B = face_dx_W_B + 0.5 * G%dyCu(I-1,j) * phi_B * P_star
      enddo

      ! ======================================================================
      ! East Face (I) of cell (i,j)
      ! Local cell is left (i), Neighbor cell is right (i+1)
      ! ======================================================================
      h_loc_A = max(CS%h_nodal(i,j,2,1), CS%min_h_shelf)
      h_loc_B = max(CS%h_nodal(i,j,2,2), CS%min_h_shelf)
      b_loc_A = bed_corners(2,1) ; b_loc_B = bed_corners(2,2)

      h_ngh_A = max(CS%h_nodal(i+1,j,1,1), CS%min_h_shelf)
      h_ngh_B = max(CS%h_nodal(i+1,j,1,2), CS%min_h_shelf)
      b_ngh_A = bed_corners(2,1) ; b_ngh_B = bed_corners(2,2)

      ! Dirichlet thickness BC: override face thicknesses with h_bdry_val on any hmask==3 side.
      loc_is_bc = (ISS%hmask(i,j) == 3)
      ngh_is_bc = (ISS%hmask(i+1,j) == 3)
      if (loc_is_bc) then
        h_loc_A = max(CS%h_bdry_val(i,j), CS%min_h_shelf) ; h_loc_B = h_loc_A
      endif
      if (ngh_is_bc) then
        h_ngh_A = max(CS%h_bdry_val(i+1,j), CS%min_h_shelf) ; h_ngh_B = h_ngh_A
      endif

      is_ext_bdry = ((CS%u_face_mask_bdry(I,j) == 2) .or. &
                    ((ISS%hmask(i+1,j) == 0 .or. ISS%hmask(i+1,j) == 2) .and. &
                     (CS%reentrant_x .or. (i+i_off /= giec))))
      is_wall = (.not. is_ext_bdry) .and. ( &
                  ((.not. CS%reentrant_x) .and. (i+i_off == giec)) .or. &
                  (int(CS%u_face_mask_bdry(I,j)) == 3) .or. &
                  (int(CS%u_face_mask_bdry(I,j)) == 5) )

      do gp_face=1,2
        t_face = xquad(gp_face)
        h_loc = (1.0 - t_face)*h_loc_A + t_face*h_loc_B
        b_loc = (1.0 - t_face)*b_loc_A + t_face*b_loc_B

        if (rhoi_rhow * h_loc - b_loc > 0.0) then
          P_loc = 0.5 * grav * rho * h_loc**2
        else
          P_loc = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_loc**2
        endif

        if (is_ext_bdry) then
          ! Ice-front Neumann (natural BC: (h*sigma_dev)*n = P_ice - P_ocean).
          ! With the flotation-branched volume IBP p_term_vol, the per-cell
          ! face contribution sums to (P_loc - P_star)*n*phi*dS, so:
          !   grounded: P_loc = 0.5*rho*g*h^2 and we set P_star = P_ocean
          !             to recover net (P_ice - P_ocean).
          !   floating: P_loc = (1 - rhoi_rhow)*0.5*rho*g*h^2 already
          !             equals (P_ice - P_ocean) for d_ocean = rhoi_rhow*h,
          !             so the IBP face P_star must vanish.
          if (rhoi_rhow * h_loc - b_loc > 0.0) then
            d_ocean = max(0.0, min(b_loc, rhoi_rhow * h_loc))
            P_star = 0.5 * grav * rhow * d_ocean**2
          else
            P_star = 0.0
          endif
        else if (is_wall) then
          ! Velocity-Dirichlet global wall: mirror p_term_vol so the IBP
          ! face term exactly cancels the local volume IBP face term. The
          ! flotation-branched P_loc equals p_term_vol per QP, so use it
          ! directly.
          P_star = P_loc
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Interior face: central P_star. See west-face block for rationale.
          if (loc_is_bc .and. .not. ngh_is_bc) then
            P_star = P_loc
          elseif (ngh_is_bc .and. .not. loc_is_bc) then
            P_star = P_ngh
          else
            P_star = 0.5 * (P_loc + P_ngh)
          endif
        endif

        phi_A = 1.0 - t_face ; phi_B = t_face
        ! Face normal n_x = 1 -> Boundary Integral = - P_star
        face_dx_E_A = face_dx_E_A - 0.5 * G%dyCu(I,j) * phi_A * P_star
        face_dx_E_B = face_dx_E_B - 0.5 * G%dyCu(I,j) * phi_B * P_star
      enddo

      ! ======================================================================
      ! South Face (J-1) of cell (i,j)
      ! Local cell is top (j), Neighbor cell is bottom (j-1)
      ! ======================================================================
      h_loc_A = max(CS%h_nodal(i,j,1,1), CS%min_h_shelf)
      h_loc_B = max(CS%h_nodal(i,j,2,1), CS%min_h_shelf)
      b_loc_A = bed_corners(1,1) ; b_loc_B = bed_corners(2,1)

      h_ngh_A = max(CS%h_nodal(i,j-1,1,2), CS%min_h_shelf)
      h_ngh_B = max(CS%h_nodal(i,j-1,2,2), CS%min_h_shelf)
      b_ngh_A = bed_corners(1,1) ; b_ngh_B = bed_corners(2,1)

      ! Dirichlet thickness BC: override face thicknesses with h_bdry_val on any hmask==3 side.
      loc_is_bc = (ISS%hmask(i,j) == 3)
      ngh_is_bc = (ISS%hmask(i,j-1) == 3)
      if (loc_is_bc) then
        h_loc_A = max(CS%h_bdry_val(i,j), CS%min_h_shelf) ; h_loc_B = h_loc_A
      endif
      if (ngh_is_bc) then
        h_ngh_A = max(CS%h_bdry_val(i,j-1), CS%min_h_shelf) ; h_ngh_B = h_ngh_A
      endif

      is_ext_bdry = ((CS%v_face_mask_bdry(i,J-1) == 2) .or. &
                    ((ISS%hmask(i,j-1) == 0 .or. ISS%hmask(i,j-1) == 2) .and. &
                     (CS%reentrant_y .or. (j+j_off /= gjsc))))
      is_wall = (.not. is_ext_bdry) .and. ( &
                  ((.not. CS%reentrant_y) .and. (j+j_off == gjsc)) .or. &
                  (int(CS%v_face_mask_bdry(i,J-1)) == 3) .or. &
                  (int(CS%v_face_mask_bdry(i,J-1)) == 5) )

      do gp_face=1,2
        t_face = xquad(gp_face)
        h_loc = (1.0 - t_face)*h_loc_A + t_face*h_loc_B
        b_loc = (1.0 - t_face)*b_loc_A + t_face*b_loc_B

        if (rhoi_rhow * h_loc - b_loc > 0.0) then
          P_loc = 0.5 * grav * rho * h_loc**2
        else
          P_loc = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_loc**2
        endif

        if (is_ext_bdry) then
          ! Ice-front Neumann (natural BC: (h*sigma_dev)*n = P_ice - P_ocean).
          ! With the flotation-branched volume IBP p_term_vol, the per-cell
          ! face contribution sums to (P_loc - P_star)*n*phi*dS, so:
          !   grounded: P_loc = 0.5*rho*g*h^2 and we set P_star = P_ocean
          !             to recover net (P_ice - P_ocean).
          !   floating: P_loc = (1 - rhoi_rhow)*0.5*rho*g*h^2 already
          !             equals (P_ice - P_ocean) for d_ocean = rhoi_rhow*h,
          !             so the IBP face P_star must vanish.
          if (rhoi_rhow * h_loc - b_loc > 0.0) then
            d_ocean = max(0.0, min(b_loc, rhoi_rhow * h_loc))
            P_star = 0.5 * grav * rhow * d_ocean**2
          else
            P_star = 0.0
          endif
        else if (is_wall) then
          ! Velocity-Dirichlet global wall: mirror p_term_vol so the IBP
          ! face term exactly cancels the local volume IBP face term. The
          ! flotation-branched P_loc equals p_term_vol per QP, so use it
          ! directly.
          P_star = P_loc
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Interior face: central P_star. See west-face block for rationale.
          if (loc_is_bc .and. .not. ngh_is_bc) then
            P_star = P_loc
          elseif (ngh_is_bc .and. .not. loc_is_bc) then
            P_star = P_ngh
          else
            P_star = 0.5 * (P_loc + P_ngh)
          endif
        endif

        phi_A = 1.0 - t_face ; phi_B = t_face
        ! Face normal n_y = -1 -> Boundary Integral = + P_star
        face_dy_S_A = face_dy_S_A + 0.5 * G%dxCv(i,J-1) * phi_A * P_star
        face_dy_S_B = face_dy_S_B + 0.5 * G%dxCv(i,J-1) * phi_B * P_star
      enddo

      ! ======================================================================
      ! North Face (J) of cell (i,j)
      ! Local cell is bottom (j), Neighbor cell is top (j+1)
      ! ======================================================================
      h_loc_A = max(CS%h_nodal(i,j,1,2), CS%min_h_shelf)
      h_loc_B = max(CS%h_nodal(i,j,2,2), CS%min_h_shelf)
      b_loc_A = bed_corners(1,2) ; b_loc_B = bed_corners(2,2)

      h_ngh_A = max(CS%h_nodal(i,j+1,1,1), CS%min_h_shelf)
      h_ngh_B = max(CS%h_nodal(i,j+1,2,1), CS%min_h_shelf)
      b_ngh_A = bed_corners(1,2) ; b_ngh_B = bed_corners(2,2)

      ! Dirichlet thickness BC: override face thicknesses with h_bdry_val on any hmask==3 side.
      loc_is_bc = (ISS%hmask(i,j) == 3)
      ngh_is_bc = (ISS%hmask(i,j+1) == 3)
      if (loc_is_bc) then
        h_loc_A = max(CS%h_bdry_val(i,j), CS%min_h_shelf) ; h_loc_B = h_loc_A
      endif
      if (ngh_is_bc) then
        h_ngh_A = max(CS%h_bdry_val(i,j+1), CS%min_h_shelf) ; h_ngh_B = h_ngh_A
      endif

      is_ext_bdry = ((CS%v_face_mask_bdry(i,J) == 2) .or. &
                    ((ISS%hmask(i,j+1) == 0 .or. ISS%hmask(i,j+1) == 2) .and. &
                     (CS%reentrant_y .or. (j+j_off /= gjec))))
      is_wall = (.not. is_ext_bdry) .and. ( &
                  ((.not. CS%reentrant_y) .and. (j+j_off == gjec)) .or. &
                  (int(CS%v_face_mask_bdry(i,J)) == 3) .or. &
                  (int(CS%v_face_mask_bdry(i,J)) == 5) )

      do gp_face=1,2
        t_face = xquad(gp_face)
        h_loc = (1.0 - t_face)*h_loc_A + t_face*h_loc_B
        b_loc = (1.0 - t_face)*b_loc_A + t_face*b_loc_B

        if (rhoi_rhow * h_loc - b_loc > 0.0) then
          P_loc = 0.5 * grav * rho * h_loc**2
        else
          P_loc = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_loc**2
        endif

        if (is_ext_bdry) then
          ! Ice-front Neumann (natural BC: (h*sigma_dev)*n = P_ice - P_ocean).
          ! With the flotation-branched volume IBP p_term_vol, the per-cell
          ! face contribution sums to (P_loc - P_star)*n*phi*dS, so:
          !   grounded: P_loc = 0.5*rho*g*h^2 and we set P_star = P_ocean
          !             to recover net (P_ice - P_ocean).
          !   floating: P_loc = (1 - rhoi_rhow)*0.5*rho*g*h^2 already
          !             equals (P_ice - P_ocean) for d_ocean = rhoi_rhow*h,
          !             so the IBP face P_star must vanish.
          if (rhoi_rhow * h_loc - b_loc > 0.0) then
            d_ocean = max(0.0, min(b_loc, rhoi_rhow * h_loc))
            P_star = 0.5 * grav * rhow * d_ocean**2
          else
            P_star = 0.0
          endif
        else if (is_wall) then
          ! Velocity-Dirichlet global wall: mirror p_term_vol so the IBP
          ! face term exactly cancels the local volume IBP face term. The
          ! flotation-branched P_loc equals p_term_vol per QP, so use it
          ! directly.
          P_star = P_loc
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Interior face: central P_star. See west-face block for rationale.
          if (loc_is_bc .and. .not. ngh_is_bc) then
            P_star = P_loc
          elseif (ngh_is_bc .and. .not. loc_is_bc) then
            P_star = P_ngh
          else
            P_star = 0.5 * (P_loc + P_ngh)
          endif
        endif

        phi_A = 1.0 - t_face ; phi_B = t_face
        ! Face normal n_y = 1 -> Boundary Integral = - P_star
        face_dy_N_A = face_dy_N_A - 0.5 * G%dxCv(i,J) * phi_A * P_star
        face_dy_N_B = face_dy_N_B - 0.5 * G%dxCv(i,J) * phi_B * P_star
      enddo

      ! Combine the per-cell volume-integral total vol_d*(m,n) with the per-face
      ! boundary contributions to get this cell's contribution to each of its 4 corner-nodes.
      cell_dx_node(1,1) = vol_dx(1,1) + face_dx_W_A  ! SW corner, W face
      cell_dy_node(1,1) = vol_dy(1,1) + face_dy_S_A  ! SW corner, S face
      cell_dx_node(2,1) = vol_dx(2,1) + face_dx_E_A  ! SE corner, E face
      cell_dy_node(2,1) = vol_dy(2,1) + face_dy_S_B  ! SE corner, S face
      cell_dx_node(1,2) = vol_dx(1,2) + face_dx_W_B  ! NW corner, W face
      cell_dy_node(1,2) = vol_dy(1,2) + face_dy_N_A  ! NW corner, N face
      cell_dx_node(2,2) = vol_dx(2,2) + face_dx_E_B  ! NE corner, E face
      cell_dy_node(2,2) = vol_dy(2,2) + face_dy_N_B  ! NE corner, N face

      taudx_b(I-1,J-1,4) = taudx_b(I-1,J-1,4) + cell_dx_node(1,1)
      taudy_b(I-1,J-1,4) = taudy_b(I-1,J-1,4) + cell_dy_node(1,1)
      taudx_b(I  ,J-1,3) = taudx_b(I  ,J-1,3) + cell_dx_node(2,1)
      taudy_b(I  ,J-1,3) = taudy_b(I  ,J-1,3) + cell_dy_node(2,1)
      taudx_b(I-1,J  ,2) = taudx_b(I-1,J  ,2) + cell_dx_node(1,2)
      taudy_b(I-1,J  ,2) = taudy_b(I-1,J  ,2) + cell_dy_node(1,2)
      taudx_b(I  ,J  ,1) = taudx_b(I  ,J  ,1) + cell_dx_node(2,2)
      taudy_b(I  ,J  ,1) = taudy_b(I  ,J  ,1) + cell_dy_node(2,2)

    endif
  enddo ; enddo

  ! Final per-node reduction: combine the 4 surrounding-cell contributions with
  ! a diagonal + off-diagonal pair-sum.
  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB
    taudx(I,J) = (taudx_b(I,J,1) + taudx_b(I,J,4)) + (taudx_b(I,J,2) + taudx_b(I,J,3))
    taudy(I,J) = (taudy_b(I,J,1) + taudy_b(I,J,4)) + (taudy_b(I,J,2) + taudy_b(I,J,3))
  enddo ; enddo

end subroutine calc_shelf_driving_stress_DG

!> Accumulate the ice-front Neumann face contribution for one face of one
!! element into the per-corner accumulators face_A, face_B. Used by the
!! strong-form driving stress at external (ocean) boundary faces.
!! Integrand at each face Gauss point:
!!   face_sign * 1/2 * face_length * phi * (P_ice - P_ocean)
!! where P_ice = 1/2 rho g h^2 and P_ocean = 1/2 rhow g d_ocean^2.
!! face_sign is the outward unit-normal component on the relevant axis
!! (-1 on W/S faces, +1 on E/N faces) and absorbs the n-dot-axis sign.
subroutine add_strong_ice_front_face(face_length, face_sign, &
    h_corner_A, h_corner_B, b_corner_A, b_corner_B, &
    h_bdry_val, loc_is_bc, rho, rhow, rhoi_rhow, grav, min_h_shelf, &
    xquad, face_A, face_B)
  real, intent(in)    :: face_length      !< Face length [L ~> m]
  real, intent(in)    :: face_sign        !< +1 or -1, outward unit-normal component [nondim]
  real, intent(in)    :: h_corner_A       !< h_nodal at face endpoint A [Z ~> m]
  real, intent(in)    :: h_corner_B       !< h_nodal at face endpoint B [Z ~> m]
  real, intent(in)    :: b_corner_A       !< bed depth at face endpoint A [Z ~> m]
  real, intent(in)    :: b_corner_B       !< bed depth at face endpoint B [Z ~> m]
  real, intent(in)    :: h_bdry_val       !< Dirichlet thickness BC value [Z ~> m]
  logical, intent(in) :: loc_is_bc        !< True to override corner h with h_bdry_val
  real, intent(in)    :: rho              !< Ice density [R ~> kg m-3]
  real, intent(in)    :: rhow             !< Ocean density [R ~> kg m-3]
  real, intent(in)    :: rhoi_rhow        !< rho / rhow [nondim]
  real, intent(in)    :: grav             !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real, intent(in)    :: min_h_shelf      !< Lower clamp on h [Z ~> m]
  real, dimension(2), intent(in) :: xquad !< 2-point Gauss-Legendre nodes on [0,1] [nondim]
  real, intent(inout) :: face_A           !< Accumulator for corner A [R L3 Z T-2 ~> kg m s-2]
  real, intent(inout) :: face_B           !< Accumulator for corner B [R L3 Z T-2 ~> kg m s-2]

  real :: h_A, h_B, b_A, b_B
  real :: t_face, h_loc, b_loc
  real :: P_ice, d_ocean, P_ocean, P_face, phi_A, phi_B
  integer :: gp_face

  h_A = max(h_corner_A, min_h_shelf)
  h_B = max(h_corner_B, min_h_shelf)
  if (loc_is_bc) then
    h_A = max(h_bdry_val, min_h_shelf) ; h_B = h_A
  endif
  b_A = b_corner_A ; b_B = b_corner_B

  do gp_face = 1, 2
    t_face = xquad(gp_face)
    h_loc = (1.0 - t_face)*h_A + t_face*h_B
    b_loc = (1.0 - t_face)*b_A + t_face*b_B
    P_ice = 0.5 * grav * rho * h_loc**2
    d_ocean = max(0.0, min(b_loc, rhoi_rhow * h_loc))
    P_ocean = 0.5 * grav * rhow * d_ocean**2
    P_face = P_ice - P_ocean
    phi_A = 1.0 - t_face ; phi_B = t_face
    face_A = face_A + face_sign * 0.5 * face_length * phi_A * P_face
    face_B = face_B + face_sign * 0.5 * face_length * phi_B * P_face
  enddo
end subroutine add_strong_ice_front_face

!> Accumulate the interior-face Dirac correction for one face of one element
!! into the per-corner accumulators face_A, face_B. Used by the strong-form
!! driving stress at interior hmask=1 / hmask=1 (or hmask=3) faces when
!! K_thresh >= 0. The broken-Q1 thickness jump at the face contributes a
!! distributional Dirac source the per-cell strong-form quadrature misses;
!! adding this term restores the correct weak-form RHS. Per face Gauss point:
!!   face_sign * 1/4 * face_length * phi * blend * rho*g*{h}*[s]
!! where {h} = 1/2*(h_loc + h_ngh), [s] = s_loc - s_ngh, and each side's
!! surface elevation s uses its own flotation test:
!!   grounded side: s = h - b
!!   floating side: s = (1 - rho_i/rho_w) * h
!! At uniformly-grounded faces this reduces to rho*g*{h}*[h] = [P]; at
!! uniformly-floating faces to (1 - rho_i/rho_w)*[P]; at mixed-flotation
!! faces it remains distributionally correct (the per-side-P difference
!! formula does not, missing terms of order 1/2*rho*g*r*h^2).
!! The 1/4 prefactor is the product of the 1/2 Gauss weight on [0,1] and
!! the 1/2 per-cell share of the inter-cell edge integral. blend = 1 when
!! K_thresh = 0 (full central edge flux); K_thresh > 0 ramps blend from 0
!! at small relative jumps to 1 at large jumps.
subroutine add_strong_mixed_interior_face(face_length, face_sign, &
    h_loc_A, h_loc_B, h_ngh_A, h_ngh_B, hg_loc_A, hg_loc_B, hg_ngh_A, hg_ngh_B, &
    b_corner_A, b_corner_B, &
    K_thresh, rho, rhow, rhoi_rhow, grav, min_h_shelf, &
    xquad, face_A, face_B)
  real, intent(in)    :: face_length      !< Face length [L ~> m]
  real, intent(in)    :: face_sign        !< +1 or -1, outward unit-normal component [nondim]
  real, intent(in)    :: h_loc_A          !< Local-cell h_nodal at face endpoint A [Z ~> m]
  real, intent(in)    :: h_loc_B          !< Local-cell h_nodal at face endpoint B [Z ~> m]
  real, intent(in)    :: h_ngh_A          !< Neighbour-cell h_nodal at face endpoint A [Z ~> m]
  real, intent(in)    :: h_ngh_B          !< Neighbour-cell h_nodal at face endpoint B [Z ~> m]
  real, intent(in)    :: hg_loc_A         !< Local-side gate thickness for the flotation test at
                                          !! endpoint A; equals h_loc_A unless DG_GL_GATE_CONTINUOUS [Z ~> m]
  real, intent(in)    :: hg_loc_B         !< Local-side gate thickness at endpoint B [Z ~> m]
  real, intent(in)    :: hg_ngh_A         !< Neighbour-side gate thickness at endpoint A; equals the
                                          !! local-side gate under DG_GL_GATE_CONTINUOUS [Z ~> m]
  real, intent(in)    :: hg_ngh_B         !< Neighbour-side gate thickness at endpoint B [Z ~> m]
  real, intent(in)    :: b_corner_A       !< bed depth at face endpoint A [Z ~> m]
  real, intent(in)    :: b_corner_B       !< bed depth at face endpoint B [Z ~> m]
  real, intent(in)    :: K_thresh         !< Venkatakrishnan-style blend threshold [nondim]
  real, intent(in)    :: rho              !< Ice density [R ~> kg m-3]
  real, intent(in)    :: rhow             !< Ocean density [R ~> kg m-3]
  real, intent(in)    :: rhoi_rhow        !< rho / rhow [nondim]
  real, intent(in)    :: grav             !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real, intent(in)    :: min_h_shelf      !< Lower clamp on h [Z ~> m]
  real, dimension(2), intent(in) :: xquad !< 2-point Gauss-Legendre nodes on [0,1] [nondim]
  real, intent(inout) :: face_A           !< Accumulator for corner A [R L3 Z T-2 ~> kg m s-2]
  real, intent(inout) :: face_B           !< Accumulator for corner B [R L3 Z T-2 ~> kg m s-2]

  real :: hL_A, hL_B, hN_A, hN_B, b_A, b_B
  real :: hgL_A, hgL_B, hgN_A, hgN_B
  real :: t_face, h_loc, h_ngh, b_loc
  real :: hg_loc, hg_ngh
  real :: s_loc, s_ngh, h_avg, jump_factor
  real :: delta_h, r2, s_blend, phi_A, phi_B
  integer :: gp_face

  hL_A = max(h_loc_A, min_h_shelf) ; hL_B = max(h_loc_B, min_h_shelf)
  hN_A = max(h_ngh_A, min_h_shelf) ; hN_B = max(h_ngh_B, min_h_shelf)
  hgL_A = max(hg_loc_A, min_h_shelf) ; hgL_B = max(hg_loc_B, min_h_shelf)
  hgN_A = max(hg_ngh_A, min_h_shelf) ; hgN_B = max(hg_ngh_B, min_h_shelf)
  b_A = b_corner_A ; b_B = b_corner_B

  do gp_face = 1, 2
    t_face = xquad(gp_face)
    h_loc = (1.0 - t_face)*hL_A + t_face*hL_B
    h_ngh = (1.0 - t_face)*hN_A + t_face*hN_B
    b_loc = (1.0 - t_face)*b_A + t_face*b_B
    hg_loc = (1.0 - t_face)*hgL_A + t_face*hgL_B
    hg_ngh = (1.0 - t_face)*hgN_A + t_face*hgN_B

    ! Per-side surface elevation. The flotation test uses the gate thickness
    ! (h_flot under DG_GL_GATE_CONTINUOUS, in which case both sides see the same
    ! gate and take the same branch); the s magnitudes use each side's own h.
    if (rhoi_rhow * hg_loc - b_loc > 0.0) then
      s_loc = h_loc - b_loc                          ! grounded
    else
      s_loc = (1.0 - rhoi_rhow) * h_loc              ! floating
    endif
    if (rhoi_rhow * hg_ngh - b_loc > 0.0) then
      s_ngh = h_ngh - b_loc                          ! grounded
    else
      s_ngh = (1.0 - rhoi_rhow) * h_ngh              ! floating
    endif

    h_avg = 0.5 * (h_loc + h_ngh)
    ! Universal Dirac integrand: rho * grav * {h} * [s]. Reduces to [P_eff]
    ! at uniformly-grounded and uniformly-floating faces; remains correct
    ! at mixed-flotation faces where [P_per-side] would miss terms.
    jump_factor = rho * grav * h_avg * (s_loc - s_ngh)

    ! Blend factor: K_thresh == 0 recovers the full distributional edge
    ! correction (blend = 1 everywhere); K_thresh > 0 ramps smoothly from
    ! strong form (small jumps) to full correction (large jumps).
    if (K_thresh == 0.0) then
      s_blend = 1.0
    else
      delta_h = h_loc - h_ngh
      r2 = (delta_h * delta_h) / max(h_avg * h_avg, 1.0e-20)
      s_blend = r2 / (r2 + K_thresh**2)
    endif

    phi_A = 1.0 - t_face ; phi_B = t_face
    face_A = face_A + face_sign * 0.25 * face_length * phi_A * s_blend * jump_factor
    face_B = face_B + face_sign * 0.25 * face_length * phi_B * s_blend * jump_factor
  enddo
end subroutine add_strong_mixed_interior_face

!> Strong-form (non-IBP) DG(1) driving stress with optional scale-aware face
!! flux at interior hmask=1 / hmask=1 faces. The volume integral is
!!   int phi * (-rho*g*h*grad(s)) dA,
!! evaluated at 2x2 Gauss points (or nsub x nsub x 2x2 sub-Gauss in the GL
!! band when GL_regularize is on) using the Q1 nodal basis for h, grad(h)
!! and corner bed_node for b, grad(b). The ice-front Neumann BC (ocean
!! back-pressure) is added as an explicit face integral
!!   int phi * (1/2*rho*g*h^2 - 1/2*rhow*g*d_ocean^2) * n dS
!! over external boundary faces. Walls contribute nothing.
!!
!! Interior face flux (optional, off by default): when
!! CS%dg_face_flux_K_thresh >= 0, a per-cell face contribution
!!   s * int phi * 0.5 * (P_loc - P_ngh) * n dS
!! is added at hmask=1/hmask=1 faces. K_thresh = 0 forces s = 1 (full
!! central-IBP face flux). K_thresh > 0 ramps via the Venkatakrishnan-style
!! smooth blend  s = r^2 / (r^2 + K^2),  r = |[h]|/h_avg, recovering central
!! flux at large r and vanishing for smooth h. K_thresh < 0 disables the
!! face flux entirely (pure strong form).
subroutine calc_shelf_driving_stress_DG_strong(CS, ISS, G, US, taudx, taudy, OD)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: taudx  !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: taudy  !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: OD  !< Ocean floor depth at tracer points [Z ~> m].

  real :: rho        ! Ice density [R ~> kg m-3]
  real :: rhow       ! Reference ocean density [R ~> kg m-3]
  real :: rhoi_rhow  ! Ice/ocean density ratio [nondim]
  real :: grav       ! Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: h_gp                  ! Ice thickness at a qp [Z ~> m]
  real :: hg_gp                 ! Gate thickness at a qp for the flotation test [Z ~> m]
  real :: dhdx_gp, dhdy_gp      ! Thickness gradients in physical coords [Z L-1 ~> nondim]
  real :: bed_gp                ! Bed depth at a qp [Z ~> m]
  real :: dbdx_gp, dbdy_gp      ! Bed-depth gradients in physical coords [Z L-1 ~> nondim]
  real :: dbdx_ref, dbdy_ref    ! Bed-depth gradients in reference coords [Z ~> m]
  real :: dsdx_gp, dsdy_gp      ! Surface gradients in physical coords [nondim]
                                ! Note: CS%bed_node is depth-positive (matches CS%bed_elev),
                                ! so for grounded ice s = h - bed and grad(s) = grad(h) - grad(bed).
  real, dimension(:,:,:,:), pointer :: hgate ! Thickness field for the flotation
                                ! test: h_flot under DG_GL_GATE_CONTINUOUS, else h_nodal [Z ~> m]
  real :: a_qp, d_qp            ! Per-qp interpolated cell-edge spacings [L ~> m]
  real :: weight                ! Per-qp quadrature weight including Jacobian [L2 ~> m2]
  real :: phi_val               ! Bilinear nodal basis value at a qp [nondim]
  real :: bed_corners(2,2)      ! Bed depth at the 4 B-grid corners of an element [Z ~> m]
  real :: dxCv_S, dxCv_N        ! Cell-edge spacings on south and north faces [L ~> m]
  real :: dyCu_W, dyCu_E        ! Cell-edge spacings on west and east faces [L ~> m]
  real :: fx_gp, fy_gp          ! Driving-force density at a qp [R Z L T-2 ~> kg m-1 s-2]

  ! Face contribution dispatch flags (per-face quadrature is delegated to
  ! add_strong_ice_front_face / add_strong_mixed_interior_face).
  logical :: is_ext_bdry         ! True if the face is an external (ocean) boundary
  logical :: loc_is_bc           ! True if local cell has hmask==3 (Dirichlet thickness BC)
  real :: face_dx_W_A, face_dx_W_B  ! West face x-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]
  real :: face_dx_E_A, face_dx_E_B  ! East face x-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]
  real :: face_dy_S_A, face_dy_S_B  ! South face y-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]
  real :: face_dy_N_A, face_dy_N_B  ! North face y-stress contrib to nodes A,B [R L3 Z T-2 ~> kg m s-2]

  real, dimension(2,2,2,2) :: qp_dx, qp_dy ! Per-QP contributions to the 4 cell-corners,
                                           ! indexed (qx,qy,m,n) [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2) :: vol_dx, vol_dy   ! Per-corner cell-volume total [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2) :: cell_dx_node, cell_dy_node ! Per-corner volume + Neumann face total
                                                     ! [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2) :: slope_x_gp, slope_y_gp ! Per-QP surface slopes pre-multiplied by
                                                 ! Jacobian a*d for physical-area weighting [Z L ~> m]
  real, dimension(2,2) :: weight_gp ! Per-QP Jacobian a*d [L2 ~> m2]
  real :: slope_w_sum  ! Sum of per-QP Jacobian weights [L2 ~> m2]

  real, dimension(SZDIB_(G),SZDJB_(G),4) :: taudx_b, taudy_b
                                           !< Per-node 4-slot driving-stress
                                           !! accumulator, slot k indexed by which
                                           !! of the 4 surrounding cells contributes
                                           !! (1=SW, 2=SE, 3=NW, 4=NE)
                                           !! [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2) :: xquad
  real :: d_corner_max, d_corner_min ! Max/min own-h flotation deficit over the 4 cell
                                     ! corners, for the subgrid dispatch test [Z ~> m]
  integer :: i, j, iq, jq, isc, iec, jsc, jec, m, n
  integer :: i_off, j_off, gisc, giec, gjsc, gjec
  logical :: calc_slope_diag
  logical :: is_grounded
  logical :: use_subgrid_cell ! True if this cell's volume integral uses subgrid quadrature

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = G%idg_offset ; j_off = G%jdg_offset
  gisc = 1 ; gjsc = 1
  giec = G%domain%niglobal ; gjec = G%domain%njglobal
  rho = CS%density_ice
  rhow = CS%density_ocean_avg
  grav = CS%g_Earth
  rhoi_rhow = rho / rhow

  xquad(1) = 0.5 * (1.0 - sqrt(1.0/3.0))
  xquad(2) = 0.5 * (1.0 + sqrt(1.0/3.0))

  taudx(:,:) = 0.0 ; taudy(:,:) = 0.0
  taudx_b(:,:,:) = 0.0 ; taudy_b(:,:,:) = 0.0

  calc_slope_diag = (CS%id_sx_shelf > 0 .or. CS%id_sy_shelf > 0 .or. CS%id_surf_slope_mag_shelf > 0)

  ! Driving-stress branch gating. By default the driving stress keeps its own-h
  ! flotation tests even under DG_GL_GATE_CONTINUOUS: per-side own-h branching is
  ! algebraically s = max(h - b, (1-r)*h), continuous in h, so the assembled force
  ! is continuous in the state. Gating the branch with h_flot instead makes s (and
  ! especially the face Dirac term rho*g*{h}*[s], which switches between [h] and
  ! (1-r)*[h]) jump discontinuously whenever the gate sign at a face flips while
  ! the local h is off flotation. That force discontinuity acts as a stiff spring
  ! pinning the steady grounding line at the configuration where the gate flips
  ! (node-average deficit = 0, i.e. GL locked at cell faces).
  hgate => CS%h_nodal
  if (CS%dg_gl_gate_continuous .and. CS%dg_gl_gate_driving_stress) hgate => CS%h_flot

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if (ISS%hmask(i,j) /= 1 .and. ISS%hmask(i,j) /= 3) cycle

    bed_corners(1,1) = CS%bed_node(I-1,J-1)
    bed_corners(2,1) = CS%bed_node(I,J-1)
    bed_corners(1,2) = CS%bed_node(I-1,J)
    bed_corners(2,2) = CS%bed_node(I,J)

    dxCv_S = G%dxCv(i,J-1) ; dxCv_N = G%dxCv(i,J)
    dyCu_W = G%dyCu(I-1,j) ; dyCu_E = G%dyCu(I,j)

    ! Volume integral: subgrid dispatch in the GL band, main-grid 2x2 Gauss otherwise.
    ! The subgrid quadrature exists to resolve the kink of THIS integrand, whose
    ! branch tests use the cell's own h_nodal; a bilinear deficit attains its
    ! extrema at the corners, so the own-h flotation contour intersects the cell
    ! iff the 4 own corner deficits have mixed signs. Under DG_GL_GATE_CONTINUOUS
    ! ground_frac is h_flot-based and can mis-dispatch (subgrid on a smooth
    ! integrand, or 2x2 Gauss on a kinked one), so the corner-sign test is used
    ! instead there. Without the gate, ground_frac is own-h-based and the frac
    ! test is kept to preserve answers.
    if (CS%dg_gl_gate_continuous) then
      d_corner_max = max( max((rhoi_rhow * max(CS%h_nodal(i,j,1,1), CS%min_h_shelf)) - bed_corners(1,1), &
                              (rhoi_rhow * max(CS%h_nodal(i,j,2,2), CS%min_h_shelf)) - bed_corners(2,2)), &
                          max((rhoi_rhow * max(CS%h_nodal(i,j,2,1), CS%min_h_shelf)) - bed_corners(2,1), &
                              (rhoi_rhow * max(CS%h_nodal(i,j,1,2), CS%min_h_shelf)) - bed_corners(1,2)) )
      d_corner_min = min( min((rhoi_rhow * max(CS%h_nodal(i,j,1,1), CS%min_h_shelf)) - bed_corners(1,1), &
                              (rhoi_rhow * max(CS%h_nodal(i,j,2,2), CS%min_h_shelf)) - bed_corners(2,2)), &
                          min((rhoi_rhow * max(CS%h_nodal(i,j,2,1), CS%min_h_shelf)) - bed_corners(2,1), &
                              (rhoi_rhow * max(CS%h_nodal(i,j,1,2), CS%min_h_shelf)) - bed_corners(1,2)) )
      use_subgrid_cell = CS%GL_regularize .and. (d_corner_max > 0.0) .and. (d_corner_min <= 0.0)
    else
      use_subgrid_cell = CS%GL_regularize .and. &
                         (CS%ground_frac(i,j) > 0.0) .and. (CS%ground_frac(i,j) < 1.0)
    endif
    if (use_subgrid_cell) then
      call calc_shelf_driving_stress_DG_strong_subgrid(CS, CS%Phisub, &
          CS%h_nodal(i,j,:,:), hgate(i,j,:,:), bed_corners, &
          dxCv_S, dxCv_N, dyCu_W, dyCu_E, &
          rho, rhow, rhoi_rhow, grav, vol_dx, vol_dy, &
          CS%sx_shelf(i,j), CS%sy_shelf(i,j), calc_slope_diag)
    else
      qp_dx(:,:,:,:) = 0.0 ; qp_dy(:,:,:,:) = 0.0
      do jq=1,2 ; do iq=1,2
        a_qp = (dxCv_S * xquad(3-jq)) + (dxCv_N * xquad(jq))
        d_qp = (dyCu_W * xquad(3-iq)) + (dyCu_E * xquad(iq))
        weight = 0.25 * (a_qp * d_qp)

        h_gp = ((CS%h_nodal(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                (CS%h_nodal(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
               ((CS%h_nodal(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                (CS%h_nodal(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
        h_gp = max(h_gp, CS%min_h_shelf)
        dhdx_gp = ( ((-xquad(3-jq))*CS%h_nodal(i,j,1,1) + ( xquad(jq))   *CS%h_nodal(i,j,2,2)) + &
                    (( xquad(3-jq))*CS%h_nodal(i,j,2,1) + (-xquad(jq))   *CS%h_nodal(i,j,1,2)) ) / a_qp
        dhdy_gp = ( ((-xquad(3-iq))*CS%h_nodal(i,j,1,1) + ( xquad(iq))   *CS%h_nodal(i,j,2,2)) + &
                    ((-xquad(iq))  *CS%h_nodal(i,j,2,1) + ( xquad(3-iq)) *CS%h_nodal(i,j,1,2)) ) / d_qp

        bed_gp = ((bed_corners(1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                  (bed_corners(2,2) * (xquad(iq)   * xquad(jq))))  + &
                 ((bed_corners(2,1) * (xquad(iq)   * xquad(3-jq))) + &
                  (bed_corners(1,2) * (xquad(3-iq) * xquad(jq))))
        dbdx_ref = ((bed_corners(1,1) * (-xquad(3-jq))) + &
                    (bed_corners(2,2) * ( xquad(jq))))  + &
                   ((bed_corners(2,1) * ( xquad(3-jq))) + &
                    (bed_corners(1,2) * (-xquad(jq))))
        dbdy_ref = ((bed_corners(1,1) * (-xquad(3-iq))) + &
                    (bed_corners(2,2) * ( xquad(iq))))  + &
                   ((bed_corners(2,1) * (-xquad(iq)))   + &
                    (bed_corners(1,2) * ( xquad(3-iq))))
        dbdx_gp = dbdx_ref / a_qp
        dbdy_gp = dbdy_ref / d_qp

        ! Flotation gate from hgate; magnitudes (h_gp, dhdx/dhdy) stay on the true
        ! DG h_nodal. hgate == h_nodal unless DG_GL_GATE_CONTINUOUS is true.
        hg_gp = ((hgate(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                 (hgate(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
                ((hgate(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                 (hgate(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
        hg_gp = max(hg_gp, CS%min_h_shelf)

        if (CS%GL_couple) then
          is_grounded = (CS%ground_frac(i,j) >= 1.0)
        else
          is_grounded = (rhoi_rhow * hg_gp - bed_gp > 0.0)
        endif

        if (is_grounded) then
          dsdx_gp = dhdx_gp - dbdx_gp
          dsdy_gp = dhdy_gp - dbdy_gp
        else
          dsdx_gp = (1.0 - rhoi_rhow) * dhdx_gp
          dsdy_gp = (1.0 - rhoi_rhow) * dhdy_gp
        endif

        fx_gp = -rho * grav * h_gp * dsdx_gp
        fy_gp = -rho * grav * h_gp * dsdy_gp

        if (calc_slope_diag) then
          weight_gp(iq,jq) = a_qp * d_qp
          slope_x_gp(iq,jq) = dsdx_gp * weight_gp(iq,jq)
          slope_y_gp(iq,jq) = dsdy_gp * weight_gp(iq,jq)
        endif

        do n=1,2 ; do m=1,2
          phi_val = (merge(xquad(iq), xquad(3-iq), m == 2)) * &
                    (merge(xquad(jq), xquad(3-jq), n == 2))
          qp_dx(iq,jq,m,n) = weight * phi_val * fx_gp
          qp_dy(iq,jq,m,n) = weight * phi_val * fy_gp
        enddo ; enddo
      enddo ; enddo

      do n=1,2 ; do m=1,2
        vol_dx(m,n) = (qp_dx(1,1,m,n) + qp_dx(2,2,m,n)) + (qp_dx(1,2,m,n) + qp_dx(2,1,m,n))
        vol_dy(m,n) = (qp_dy(1,1,m,n) + qp_dy(2,2,m,n)) + (qp_dy(1,2,m,n) + qp_dy(2,1,m,n))
      enddo ; enddo

      if (calc_slope_diag) then
        slope_w_sum = (weight_gp(1,1)+weight_gp(2,2)) + (weight_gp(1,2)+weight_gp(2,1))
        CS%sx_shelf(i,j) = ((slope_x_gp(1,1)+slope_x_gp(2,2)) + (slope_x_gp(1,2)+slope_x_gp(2,1))) / slope_w_sum
        CS%sy_shelf(i,j) = ((slope_y_gp(1,1)+slope_y_gp(2,2)) + (slope_y_gp(1,2)+slope_y_gp(2,1))) / slope_w_sum
      endif
    endif

    ! Face contributions. Two cases per face:
    ! (1) is_ext_bdry: ice-front Neumann
    !       +integral phi * (P_ice - P_ocean) * n dS
    !     with the natural n sign (-W, +E, -S, +N).
    ! (2) Interior hmask=1/hmask=1 face and CS%dg_face_flux_K_thresh >= 0:
    !     mixed-form scale-aware face flux
    !       s * integral phi * 0.5*(P_loc - P_ngh) * n dS
    !     where s = 1 when K_thresh = 0 (full central-IBP face flux), and
    !     s = r^2/(r^2 + K^2), r = |[h]|/h_avg, when K_thresh > 0 (smooth
    !     ramp from strong to central). The per-cell factor 0.5 is the
    !     leftover after the IBP-reverse identity is applied to the strong
    !     form (each adjacent cell contributes its own 0.5*(P_K - P_K')*n;
    !     the K and K' contributions sum to int phi*[P]*n at the shared
    !     B-node). K_thresh < 0 (default) disables the face flux entirely
    !     and recovers the pure strong form.
    ! Walls and other (hmask=3 / hmask=0 interior) faces contribute nothing.
    face_dx_W_A = 0.0 ; face_dx_W_B = 0.0
    face_dx_E_A = 0.0 ; face_dx_E_B = 0.0
    face_dy_S_A = 0.0 ; face_dy_S_B = 0.0
    face_dy_N_A = 0.0 ; face_dy_N_B = 0.0
    loc_is_bc = (ISS%hmask(i,j) == 3)

    ! West face (n_x = -1).
    is_ext_bdry = ((CS%u_face_mask_bdry(I-1,j) == 2) .or. &
                  ((ISS%hmask(i-1,j) == 0 .or. ISS%hmask(i-1,j) == 2) .and. &
                   (CS%reentrant_x .or. (i+i_off /= gisc))))
    if (is_ext_bdry) then
      call add_strong_ice_front_face(G%dyCu(I-1,j), -1.0, &
        CS%h_nodal(i,j,1,1), CS%h_nodal(i,j,1,2), bed_corners(1,1), bed_corners(1,2), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dx_W_A, face_dx_W_B)
    elseif (CS%dg_face_flux_K_thresh >= 0.0 .and. ISS%hmask(i-1,j) == 1.0) then
      call add_strong_mixed_interior_face(G%dyCu(I-1,j), -1.0, &
        CS%h_nodal(i,j,1,1),   CS%h_nodal(i,j,1,2), &
        CS%h_nodal(i-1,j,2,1), CS%h_nodal(i-1,j,2,2), &
        hgate(i,j,1,1),   hgate(i,j,1,2), &
        hgate(i-1,j,2,1), hgate(i-1,j,2,2), &
        bed_corners(1,1), bed_corners(1,2), CS%dg_face_flux_K_thresh, &
        rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dx_W_A, face_dx_W_B)
    endif

    ! East face (n_x = +1).
    is_ext_bdry = ((CS%u_face_mask_bdry(I,j) == 2) .or. &
                  ((ISS%hmask(i+1,j) == 0 .or. ISS%hmask(i+1,j) == 2) .and. &
                   (CS%reentrant_x .or. (i+i_off /= giec))))
    if (is_ext_bdry) then
      call add_strong_ice_front_face(G%dyCu(I,j), +1.0, &
        CS%h_nodal(i,j,2,1), CS%h_nodal(i,j,2,2), bed_corners(2,1), bed_corners(2,2), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dx_E_A, face_dx_E_B)
    elseif (CS%dg_face_flux_K_thresh >= 0.0 .and. ISS%hmask(i+1,j) == 1.0) then
      call add_strong_mixed_interior_face(G%dyCu(I,j), +1.0, &
        CS%h_nodal(i,j,2,1),   CS%h_nodal(i,j,2,2), &
        CS%h_nodal(i+1,j,1,1), CS%h_nodal(i+1,j,1,2), &
        hgate(i,j,2,1),   hgate(i,j,2,2), &
        hgate(i+1,j,1,1), hgate(i+1,j,1,2), &
        bed_corners(2,1), bed_corners(2,2), CS%dg_face_flux_K_thresh, &
        rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dx_E_A, face_dx_E_B)
    endif

    ! South face (n_y = -1).
    is_ext_bdry = ((CS%v_face_mask_bdry(i,J-1) == 2) .or. &
                  ((ISS%hmask(i,j-1) == 0 .or. ISS%hmask(i,j-1) == 2) .and. &
                   (CS%reentrant_y .or. (j+j_off /= gjsc))))
    if (is_ext_bdry) then
      call add_strong_ice_front_face(G%dxCv(i,J-1), -1.0, &
        CS%h_nodal(i,j,1,1), CS%h_nodal(i,j,2,1), bed_corners(1,1), bed_corners(2,1), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dy_S_A, face_dy_S_B)
    elseif (CS%dg_face_flux_K_thresh >= 0.0 .and. ISS%hmask(i,j-1) == 1.0) then
      call add_strong_mixed_interior_face(G%dxCv(i,J-1), -1.0, &
        CS%h_nodal(i,j,1,1),   CS%h_nodal(i,j,2,1), &
        CS%h_nodal(i,j-1,1,2), CS%h_nodal(i,j-1,2,2), &
        hgate(i,j,1,1),   hgate(i,j,2,1), &
        hgate(i,j-1,1,2), hgate(i,j-1,2,2), &
        bed_corners(1,1), bed_corners(2,1), CS%dg_face_flux_K_thresh, &
        rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dy_S_A, face_dy_S_B)
    endif

    ! North face (n_y = +1).
    is_ext_bdry = ((CS%v_face_mask_bdry(i,J) == 2) .or. &
                  ((ISS%hmask(i,j+1) == 0 .or. ISS%hmask(i,j+1) == 2) .and. &
                   (CS%reentrant_y .or. (j+j_off /= gjec))))
    if (is_ext_bdry) then
      call add_strong_ice_front_face(G%dxCv(i,J), +1.0, &
        CS%h_nodal(i,j,1,2), CS%h_nodal(i,j,2,2), bed_corners(1,2), bed_corners(2,2), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dy_N_A, face_dy_N_B)
    elseif (CS%dg_face_flux_K_thresh >= 0.0 .and. ISS%hmask(i,j+1) == 1.0) then
      call add_strong_mixed_interior_face(G%dxCv(i,J), +1.0, &
        CS%h_nodal(i,j,1,2),   CS%h_nodal(i,j,2,2), &
        CS%h_nodal(i,j+1,1,1), CS%h_nodal(i,j+1,2,1), &
        hgate(i,j,1,2),   hgate(i,j,2,2), &
        hgate(i,j+1,1,1), hgate(i,j+1,2,1), &
        bed_corners(1,2), bed_corners(2,2), CS%dg_face_flux_K_thresh, &
        rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dy_N_A, face_dy_N_B)
    endif

    ! Combine cell-volume integral with face contributions (ice-front Neumann
    ! and, when enabled, mixed-form scale-aware interior face flux).
    cell_dx_node(1,1) = vol_dx(1,1) + face_dx_W_A
    cell_dy_node(1,1) = vol_dy(1,1) + face_dy_S_A
    cell_dx_node(2,1) = vol_dx(2,1) + face_dx_E_A
    cell_dy_node(2,1) = vol_dy(2,1) + face_dy_S_B
    cell_dx_node(1,2) = vol_dx(1,2) + face_dx_W_B
    cell_dy_node(1,2) = vol_dy(1,2) + face_dy_N_A
    cell_dx_node(2,2) = vol_dx(2,2) + face_dx_E_B
    cell_dy_node(2,2) = vol_dy(2,2) + face_dy_N_B

    taudx_b(I-1,J-1,4) = taudx_b(I-1,J-1,4) + cell_dx_node(1,1)
    taudy_b(I-1,J-1,4) = taudy_b(I-1,J-1,4) + cell_dy_node(1,1)
    taudx_b(I  ,J-1,3) = taudx_b(I  ,J-1,3) + cell_dx_node(2,1)
    taudy_b(I  ,J-1,3) = taudy_b(I  ,J-1,3) + cell_dy_node(2,1)
    taudx_b(I-1,J  ,2) = taudx_b(I-1,J  ,2) + cell_dx_node(1,2)
    taudy_b(I-1,J  ,2) = taudy_b(I-1,J  ,2) + cell_dy_node(1,2)
    taudx_b(I  ,J  ,1) = taudx_b(I  ,J  ,1) + cell_dx_node(2,2)
    taudy_b(I  ,J  ,1) = taudy_b(I  ,J  ,1) + cell_dy_node(2,2)
  enddo ; enddo

  do J=G%JsdB,G%JedB ; do I=G%IsdB,G%IedB
    taudx(I,J) = (taudx_b(I,J,1) + taudx_b(I,J,4)) + (taudx_b(I,J,2) + taudx_b(I,J,3))
    taudy(I,J) = (taudy_b(I,J,1) + taudy_b(I,J,4)) + (taudy_b(I,J,2) + taudy_b(I,J,3))
  enddo ; enddo

end subroutine calc_shelf_driving_stress_DG_strong

!> Strong-form (non-IBP) subgrid volume integral of the driving stress for
!! the DG path, used in the grounding-line band when GL_regularize is on.
!! Evaluates int phi * (-rho*g*h*grad(s)) dA over nsub x nsub sub-cells,
!! with per-sub-QP floating/grounded selector via the local h_gp vs bed_gp
!! comparison. No face integrals here; the ice-front Neumann is applied at
!! the main-grid cell edges by the caller.
subroutine calc_shelf_driving_stress_DG_strong_subgrid(CS, Phisub, &
    h_nodal_cell, h_gate_cell, bed_corners, &
    dxCv_S, dxCv_N, dyCu_W, dyCu_E, &
    rho, rhow, rhoi_rhow, grav, vol_dx, vol_dy, sx_shelf, sy_shelf, calc_slope_diag)
  type(ice_shelf_dyn_CS), intent(in) :: CS    !< Ice shelf control structure
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Sub-grid quadrature weights [nondim]
  real, dimension(2,2), intent(in) :: h_nodal_cell !< Q1 nodal thickness at the 4 corners [Z ~> m]
  real, dimension(2,2), intent(in) :: h_gate_cell !< Gate thickness for the flotation test at the 4
                                          !! corners; equals h_nodal_cell unless
                                          !! DG_GL_GATE_CONTINUOUS is true [Z ~> m]
  real, dimension(2,2), intent(in) :: bed_corners !< Bed depth at the 4 cell corners [Z ~> m]
  real, intent(in) :: dxCv_S         !< Cell x-length on south face [L ~> m]
  real, intent(in) :: dxCv_N         !< Cell x-length on north face [L ~> m]
  real, intent(in) :: dyCu_W         !< Cell y-length on west face [L ~> m]
  real, intent(in) :: dyCu_E         !< Cell y-length on east face [L ~> m]
  real, intent(in) :: rho            !< Ice density [R ~> kg m-3]
  real, intent(in) :: rhow           !< Ocean density [R ~> kg m-3]
  real, intent(in) :: rhoi_rhow      !< rho/rhow [nondim]
  real, intent(in) :: grav           !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real, dimension(2,2), intent(out) :: vol_dx !< Per-corner x volume integral [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2), intent(out) :: vol_dy !< Per-corner y volume integral [R L3 Z T-2 ~> kg m s-2]
  real, intent(inout) :: sx_shelf !< Cell-average x surface slope [Z L-1 ~> nondim]
  real, intent(inout) :: sy_shelf !< Cell-average y surface slope [Z L-1 ~> nondim]
  logical :: calc_slope_diag !< True if slope diagnostics will be calculated

  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: contr_sub_dx, contr_sub_dy
  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: slope_x_gp, slope_y_gp
  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: weight_gp
  real, dimension(2,2) :: slope_x, slope_y
  real, dimension(2,2) :: weight_sum_qp
  real :: slope_w_total
  real, dimension(2,2,2,2) :: qp_dx, qp_dy
  real :: h_gp, dhdx_gp, dhdy_gp
  real :: hg_gp  ! Gate thickness at the sub-qp for the flotation test [Z ~> m]
  real :: bed_gp, dbdx_gp, dbdy_gp
  real :: dbdx_ref, dbdy_ref
  real :: dsdx_gp, dsdy_gp
  real :: fx_gp, fy_gp
  real :: y_marginal_1, y_marginal_2, x_marginal_1, x_marginal_2
  real :: a, d, weight, subarea
  integer :: nsub, i, j, qx, qy, m, n
  logical :: is_grounded

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  do j=1,nsub ; do i=1,nsub
    qp_dx(:,:,:,:) = 0.0 ; qp_dy(:,:,:,:) = 0.0
    do qy=1,2 ; do qx=1,2

      y_marginal_1 = Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1)
      y_marginal_2 = Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)
      x_marginal_1 = Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2)
      x_marginal_2 = Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)

      h_gp = ((Phisub(qx,qy,i,j,1,1)*h_nodal_cell(1,1)) + (Phisub(qx,qy,i,j,2,2)*h_nodal_cell(2,2))) + &
             ((Phisub(qx,qy,i,j,1,2)*h_nodal_cell(1,2)) + (Phisub(qx,qy,i,j,2,1)*h_nodal_cell(2,1)))
      h_gp = max(h_gp, CS%min_h_shelf)
      bed_gp = ((Phisub(qx,qy,i,j,1,1)*bed_corners(1,1)) + (Phisub(qx,qy,i,j,2,2)*bed_corners(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*bed_corners(1,2)) + (Phisub(qx,qy,i,j,2,1)*bed_corners(2,1)))

      a = (dxCv_S * y_marginal_1) + (dxCv_N * y_marginal_2)
      d = (dyCu_W * x_marginal_1) + (dyCu_E * x_marginal_2)
      weight = 0.25 * subarea * (a * d)

      dhdx_gp = ( ((-y_marginal_1) * h_nodal_cell(1,1) + ( y_marginal_2) * h_nodal_cell(2,2)) + &
                  (( y_marginal_1) * h_nodal_cell(2,1) + (-y_marginal_2) * h_nodal_cell(1,2)) ) / a
      dhdy_gp = ( ((-x_marginal_1) * h_nodal_cell(1,1) + ( x_marginal_2) * h_nodal_cell(2,2)) + &
                  ((-x_marginal_2) * h_nodal_cell(2,1) + ( x_marginal_1) * h_nodal_cell(1,2)) ) / d

      dbdx_ref = ((bed_corners(1,1) * (-y_marginal_1))  + &
                  (bed_corners(2,2) * ( y_marginal_2))) + &
                 ((bed_corners(2,1) * ( y_marginal_1))  + &
                  (bed_corners(1,2) * (-y_marginal_2)))
      dbdy_ref = ((bed_corners(1,1) * (-x_marginal_1))  + &
                  (bed_corners(2,2) * ( x_marginal_2))) + &
                 ((bed_corners(2,1) * (-x_marginal_2))  + &
                  (bed_corners(1,2) * ( x_marginal_1)))
      dbdx_gp = dbdx_ref / a
      dbdy_gp = dbdy_ref / d

      ! Sub-QP-local floating/grounded selector. For the strong form we use
      ! the local hydrostatic test even when CS%GL_couple is true, because
      ! the whole point of the subgrid is to resolve sub-cell GL position
      ! that a cell-mean ground_frac smears. The test uses the gate thickness
      ! (h_flot under DG_GL_GATE_CONTINUOUS); magnitudes keep h_nodal_cell.
      hg_gp = ((Phisub(qx,qy,i,j,1,1)*h_gate_cell(1,1)) + (Phisub(qx,qy,i,j,2,2)*h_gate_cell(2,2))) + &
              ((Phisub(qx,qy,i,j,1,2)*h_gate_cell(1,2)) + (Phisub(qx,qy,i,j,2,1)*h_gate_cell(2,1)))
      hg_gp = max(hg_gp, CS%min_h_shelf)
      is_grounded = (rhoi_rhow * hg_gp - bed_gp > 0.0)

      if (is_grounded) then
        dsdx_gp = dhdx_gp - dbdx_gp
        dsdy_gp = dhdy_gp - dbdy_gp
      else
        dsdx_gp = (1.0 - rhoi_rhow) * dhdx_gp
        dsdy_gp = (1.0 - rhoi_rhow) * dhdy_gp
      endif

      fx_gp = -rho * grav * h_gp * dsdx_gp
      fy_gp = -rho * grav * h_gp * dsdy_gp

      if (calc_slope_diag) then
        weight_gp(i,j,qx,qy) = a * d
        slope_x_gp(i,j,qx,qy) = dsdx_gp * weight_gp(i,j,qx,qy)
        slope_y_gp(i,j,qx,qy) = dsdy_gp * weight_gp(i,j,qx,qy)
      endif

      do n=1,2 ; do m=1,2
        qp_dx(qx,qy,m,n) = weight * Phisub(qx,qy,i,j,m,n) * fx_gp
        qp_dy(qx,qy,m,n) = weight * Phisub(qx,qy,i,j,m,n) * fy_gp
      enddo ; enddo
    enddo ; enddo

    do n=1,2 ; do m=1,2
      contr_sub_dx(i,j,m,n) = (qp_dx(1,1,m,n) + qp_dx(2,2,m,n)) + (qp_dx(1,2,m,n) + qp_dx(2,1,m,n))
      contr_sub_dy(i,j,m,n) = (qp_dy(1,1,m,n) + qp_dy(2,2,m,n)) + (qp_dy(1,2,m,n) + qp_dy(2,1,m,n))
    enddo ; enddo
  enddo ; enddo

  do n=1,2 ; do m=1,2
    call sum_square_matrix(vol_dx(m,n), contr_sub_dx(:,:,m,n), nsub)
    call sum_square_matrix(vol_dy(m,n), contr_sub_dy(:,:,m,n), nsub)
  enddo ; enddo

  if (calc_slope_diag) then
    do qy=1,2 ; do qx=1,2
      call sum_square_matrix(slope_x(qx,qy), slope_x_gp(:,:,qx,qy), nsub)
      call sum_square_matrix(slope_y(qx,qy), slope_y_gp(:,:,qx,qy), nsub)
      call sum_square_matrix(weight_sum_qp(qx,qy), weight_gp(:,:,qx,qy), nsub)
    enddo ; enddo
    slope_w_total = (weight_sum_qp(1,1)+weight_sum_qp(2,2)) + (weight_sum_qp(1,2)+weight_sum_qp(2,1))
    sx_shelf = ((slope_x(1,1)+slope_x(2,2)) + (slope_x(1,2)+slope_x(2,1))) / slope_w_total
    sy_shelf = ((slope_y(1,1)+slope_y(2,2)) + (slope_y(1,2)+slope_y(2,1))) / slope_w_total
  endif

end subroutine calc_shelf_driving_stress_DG_strong_subgrid

!> Subgrid GL-band volume integral of the driving stress for the DG path.
!! Evaluates the unified integration-by-parts weak form over nsub x nsub sub-cells.
subroutine calc_shelf_driving_stress_DG_subgrid(CS, Phisub, &
    h_shelf_cell, h_nodal_cell, bed_corners, &
    dxCv_S, dxCv_N, dyCu_W, dyCu_E, &
    rho, rhow, rhoi_rhow, grav, vol_dx, vol_dy, sx_shelf, sy_shelf, calc_slope_diag)
  type(ice_shelf_dyn_CS), intent(in) :: CS    !< Ice shelf control structure
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Sub-grid quadrature weights [nondim]
  real, intent(in) :: h_shelf_cell   !< Cell-averaged ice thickness [Z ~> m]
  real, dimension(2,2), intent(in) :: h_nodal_cell !< Q1 nodal thickness at the 4 corners [Z ~> m]
  real, dimension(2,2), intent(in) :: bed_corners !< Bed elevation at the 4 cell corners [Z ~> m]
  real, intent(in) :: dxCv_S         !< Cell x-length on south face [L ~> m]
  real, intent(in) :: dxCv_N         !< Cell x-length on north face [L ~> m]
  real, intent(in) :: dyCu_W         !< Cell y-length on west face [L ~> m]
  real, intent(in) :: dyCu_E         !< Cell y-length on east face [L ~> m]
  real, intent(in) :: rho            !< Ice density [R ~> kg m-3]
  real, intent(in) :: rhow           !< Ocean density [R ~> kg m-3]
  real, intent(in) :: rhoi_rhow      !< rho/rhow [nondim]
  real, intent(in) :: grav           !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real, dimension(2,2), intent(out) :: vol_dx !< Per-corner x volume integral [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2), intent(out) :: vol_dy !< Per-corner y volume integral [R L3 Z T-2 ~> kg m s-2]
  real, intent(inout) :: sx_shelf !< The cell-average x surface slope [Z L-1 ~> nondim]
  real, intent(inout) :: sy_shelf !< The cell-average y surface slope [Z L-1 ~> nondim]
  logical :: calc_slope_diag !< True if slope diagnostics will be calculated

  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: contr_sub_dx, contr_sub_dy
  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: slope_x_gp, slope_y_gp ! Per-QP surf slopes
                                                  ! pre-multiplied by Jacobian a*d for physical-area
                                                  ! weighting [Z L ~> m]
  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: weight_gp ! Per-QP Jacobian a*d used for
                                                  ! area-weighting slopes [L2 ~> m2]
  real, dimension(2,2) :: slope_x, slope_y ! slope sums for qps with same position within subcells [Z L ~> m]
  real, dimension(2,2) :: weight_sum_qp ! Jacobian weight sums per (qx,qy) [L2 ~> m2]
  real :: slope_w_total ! Total Jacobian weight summed over all sub-QPs [L2 ~> m2]
  real, dimension(2,2,2,2) :: qp_dx, qp_dy
  real :: xi_sub, eta_sub   ! DG reference coords at sub-qp ([-0.5,0.5]) [nondim]
  real :: h_gp              ! Ice thickness at sub-qp [Z ~> m]
  real :: dhdx_gp, dhdy_gp  ! Thickness gradients at sub-qp, physical coords [Z L-1 ~> nondim]
  real :: bed_gp            ! Bed elevation at sub-qp [Z ~> m]
  real :: dbdx_gp, dbdy_gp  ! Bed gradients at sub-qp, physical coords [Z L-1 ~> nondim]
  real :: dbdx_ref, dbdy_ref ! Bed gradients in reference coords (pre-Jacobian) [Z ~> m]
  real :: bottom_force_x, bottom_force_y ! Bed/water bottom drag forces [R Z L2 T-2 ~> kg s-2]
  real :: dphi_dx_ref, dphi_dy_ref ! Basis function derivatives in reference coordinates [nondim]
  real :: dphi_dx, dphi_dy  ! Basis function derivatives in physical coordinates [L-1 ~> m-1]
  real :: y_marginal_1, y_marginal_2, x_marginal_1, x_marginal_2 ! Marginal sums [nondim]
  real :: p_term_vol        ! Integrated-by-parts volume pressure term [R Z L2 T-2 ~> kg s-2]
  real :: a, d              ! Per-sub-qp interpolated cell-edge spacings [L ~> m]
  real :: weight            ! Per-sub-qp quadrature weight [L2 ~> m2]
  real :: subarea           ! 1/nsub^2 [nondim]
  integer :: nsub, i, j, qx, qy, m, n

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  do j=1,nsub ; do i=1,nsub
    qp_dx(:,:,:,:) = 0.0 ; qp_dy(:,:,:,:) = 0.0
    do qy=1,2 ; do qx=1,2

      y_marginal_1 = Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1)
      y_marginal_2 = Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)
      x_marginal_1 = Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2)
      x_marginal_2 = Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)

      ! Reference coords at the sub-qp: xi_sub = a_right(qx,i) - 0.5; marginal
      ! sum of Phisub over k for l=2 gives a_right(qx,i).
      xi_sub  = x_marginal_2 - 0.5
      eta_sub = y_marginal_2 - 0.5

      ! Nodal Q1 evaluation of h at the sub-QP using the Phisub corner-basis
      ! weights (rotation-paired). Equivalent to bed_gp's contraction.
      h_gp = ((Phisub(qx,qy,i,j,1,1)*h_nodal_cell(1,1)) + (Phisub(qx,qy,i,j,2,2)*h_nodal_cell(2,2))) + &
             ((Phisub(qx,qy,i,j,1,2)*h_nodal_cell(1,2)) + (Phisub(qx,qy,i,j,2,1)*h_nodal_cell(2,1)))
      h_gp = max(h_gp, CS%min_h_shelf)

      ! Bed at sub-qp: rotation-paired Phisub contraction of bed_corners.
      bed_gp = ((Phisub(qx,qy,i,j,1,1)*bed_corners(1,1)) + (Phisub(qx,qy,i,j,2,2)*bed_corners(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*bed_corners(1,2)) + (Phisub(qx,qy,i,j,2,1)*bed_corners(2,1)))

      ! Per-sub-qp metric via Phisub marginal sums (same pattern as CG_action_subgrid_basal).
      a = (dxCv_S * y_marginal_1) + (dxCv_N * y_marginal_2)
      d = (dyCu_W * x_marginal_1) + (dyCu_E * x_marginal_2)
      weight = 0.25 * subarea * (a * d)

      ! Nodal Q1 gradients at sub-qp via Phisub marginals; same formula structure
      ! as bed gradients below.
      dhdx_gp = ( ((-y_marginal_1) * h_nodal_cell(1,1) + ( y_marginal_2) * h_nodal_cell(2,2)) + &
                  (( y_marginal_1) * h_nodal_cell(2,1) + (-y_marginal_2) * h_nodal_cell(1,2)) ) / a
      dhdy_gp = ( ((-x_marginal_1) * h_nodal_cell(1,1) + ( x_marginal_2) * h_nodal_cell(2,2)) + &
                  ((-x_marginal_2) * h_nodal_cell(2,1) + ( x_marginal_1) * h_nodal_cell(1,2)) ) / d

      ! Reference-coord bed gradients: derivative of bilinear corner-basis at sub-qp.
      dbdx_ref = ((bed_corners(1,1) * (-y_marginal_1))  + &
                  (bed_corners(2,2) * ( y_marginal_2))) + &
                 ((bed_corners(2,1) * ( y_marginal_1))  + &
                  (bed_corners(1,2) * (-y_marginal_2)))
      dbdy_ref = ((bed_corners(1,1) * (-x_marginal_1))  + &
                  (bed_corners(2,2) * ( x_marginal_2))) + &
                 ((bed_corners(2,1) * (-x_marginal_2))  + &
                  (bed_corners(1,2) * ( x_marginal_1)))
      dbdx_gp = dbdx_ref / a
      dbdy_gp = dbdy_ref / d

      ! Flotation-branched IBP decomposition (see main-grid path for the
      ! derivation). grounded: P = 0.5*rho*g*h^2, bottom_force = rho*g*h*grad(b).
      ! floating: P = (1 - rhoi_rhow)*0.5*rho*g*h^2, bottom_force = 0.
      if (rhoi_rhow * h_gp - bed_gp > 0.0) then
        p_term_vol = 0.5 * rho * grav * h_gp**2
        bottom_force_x = rho * grav * h_gp * dbdx_gp
        bottom_force_y = rho * grav * h_gp * dbdy_gp
      else
        p_term_vol = 0.5 * (1.0 - rhoi_rhow) * rho * grav * h_gp**2
        bottom_force_x = 0.0
        bottom_force_y = 0.0
      endif

      ! For slope diagnostics. Pre-multiply by the per-sub-QP Jacobian a*d so the
      ! sub-cell sum is a physical-area integral; divide by sum of Jacobians at
      ! the end to get the area-weighted mean.
      if (calc_slope_diag) then
        weight_gp(i,j,qx,qy) = a * d
        if (rhoi_rhow * h_gp - bed_gp <= 0.0) then
          ! Floating: bottom force is water pressure on the sloped draft
          slope_x_gp(i,j,qx,qy) = (1.0 - rhoi_rhow) * dhdx_gp * weight_gp(i,j,qx,qy)
          slope_y_gp(i,j,qx,qy) = (1.0 - rhoi_rhow) * dhdy_gp * weight_gp(i,j,qx,qy)
        else
          slope_x_gp(i,j,qx,qy) = (dhdx_gp - dbdx_gp) * weight_gp(i,j,qx,qy)
          slope_y_gp(i,j,qx,qy) = (dhdy_gp - dbdy_gp) * weight_gp(i,j,qx,qy)
        endif
      endif

      ! Unified Weak-form Volume Integration applied to subgrid; p_term_vol
      ! set above by the flotation branch.
      do n=1,2 ; do m=1,2
        dphi_dx_ref = merge(1.0, -1.0, m==2) * merge(y_marginal_2, y_marginal_1, n==2)
        dphi_dy_ref = merge(x_marginal_2, x_marginal_1, m==2) * merge(1.0, -1.0, n==2)
        dphi_dx = dphi_dx_ref / a
        dphi_dy = dphi_dy_ref / d

        ! Geometric correction for the IBP pressure term on non-rectangular (e.g. lat/lon)
        ! elements. On a bilinear element, a(eta) = dxCv_S*(1-eta) + dxCv_N*eta varies with
        ! eta, so the reference-space IBP of P*d(phi)/deta requires an extra term
        ! P*phi*da/deta to satisfy the discrete divergence theorem. Without it, constant P
        ! gives a spurious metric-dependent driving stress on lat/lon grids. The 0.25*subarea
        ! factor equals weight/(a*d), converting the reference-space integral to the nodal sum.
        qp_dx(qx,qy,m,n) = (weight * dphi_dx * p_term_vol) &
                          + (weight * Phisub(qx,qy,i,j,m,n) * bottom_force_x) &
                          + (0.25 * subarea * Phisub(qx,qy,i,j,m,n) * p_term_vol * (dyCu_E - dyCu_W))
        qp_dy(qx,qy,m,n) = (weight * dphi_dy * p_term_vol) &
                          + (weight * Phisub(qx,qy,i,j,m,n) * bottom_force_y) &
                          + (0.25 * subarea * Phisub(qx,qy,i,j,m,n) * p_term_vol * (dxCv_N - dxCv_S))
      enddo ; enddo

    enddo ; enddo

    do n=1,2 ; do m=1,2
      contr_sub_dx(i,j,m,n) = (qp_dx(1,1,m,n) + qp_dx(2,2,m,n)) + &
                              (qp_dx(1,2,m,n) + qp_dx(2,1,m,n))
      contr_sub_dy(i,j,m,n) = (qp_dy(1,1,m,n) + qp_dy(2,2,m,n)) + &
                              (qp_dy(1,2,m,n) + qp_dy(2,1,m,n))
    enddo ; enddo
  enddo ; enddo

  do n=1,2 ; do m=1,2
    call sum_square_matrix(vol_dx(m,n), contr_sub_dx(:,:,m,n), nsub)
    call sum_square_matrix(vol_dy(m,n), contr_sub_dy(:,:,m,n), nsub)
  enddo ; enddo

  if (calc_slope_diag) then
    do qy=1,2 ; do qx=1,2
      call sum_square_matrix(slope_x(qx,qy), slope_x_gp(:,:,qx,qy), nsub)
      call sum_square_matrix(slope_y(qx,qy), slope_y_gp(:,:,qx,qy), nsub)
      call sum_square_matrix(weight_sum_qp(qx,qy), weight_gp(:,:,qx,qy), nsub)
    enddo ; enddo

    ! Physical-area-weighted mean: slope_*_gp already premultiplied by a*d.
    ! The 0.25 * subarea factor cancels between numerator and denominator.
    slope_w_total = (weight_sum_qp(1,1)+weight_sum_qp(2,2)) + (weight_sum_qp(1,2)+weight_sum_qp(2,1))
    sx_shelf = ((slope_x(1,1)+slope_x(2,2)) + (slope_x(1,2)+slope_x(2,1))) / slope_w_total
    sy_shelf = ((slope_y(1,1)+slope_y(2,2)) + (slope_y(1,2)+slope_y(2,1))) / slope_w_total
  endif

end subroutine calc_shelf_driving_stress_DG_subgrid

!> Reconstruct node-based bed elevation from cell-averaged bed_elev using a
!! matrix-free iterative method. The resulting bed_node values are such that
!! the bilinear interpolant over each cell integrates exactly to the cell-averaged
!! bed_elev, and the representation is continuous across cell boundaries.
!! This should be called once at initialization.
subroutine reconstruct_bed_to_nodes(CS, G, hmask)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf

  ! Local variables
  real, dimension(SZDIB_(G),SZDJB_(G)) :: bed_node_new  ! Updated node values [Z ~> m]
  real    :: residual     ! Cell constraint residual [Z ~> m]
  real    :: bilinear_avg ! Average of 4 corner values [Z ~> m]
  real    :: max_err      ! Maximum residual across all cells [Z ~> m]
  real    :: tol          ! Convergence tolerance [Z ~> m]
  real    :: c00, c10, c01, c11 ! Per-cell-of-node contributions, used to assemble
                                ! a rotation-invariant 4-element reduction [Z ~> m]
  real :: damping ! Damps the step size each iteration [nondim]

  integer :: i, j, iter, num_cells
  integer :: isc, iec, jsc, jec, IsdB, IedB, JsdB, JedB
  integer, parameter :: max_iter = 1000
  character(len=160) :: mesg  ! The text of a MOM message

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  tol = 1.0e-12

  ! Step 1: Initialize bed_node as average of surrounding cell-averaged bed_elev values
  ! (same approach as interpolate_H_to_B). The 4 cell contributions are summed with
  ! a diagonal + off-diagonal reduction so the result is invariant under 90 deg
  ! grid rotation, which cyclically permutes the four neighbour cells.
  CS%bed_node(:,:) = 0.0
  do J=jsc-1,jec ; do I=isc-1,iec
    num_cells = 0
    c00 = 0.0 ; c10 = 0.0 ; c01 = 0.0 ; c11 = 0.0
    if (hmask(I,  J  ) == 1.0 .or. hmask(I,  J  ) == 3.0) then
      c00 = CS%bed_elev(I,  J  ) ; num_cells = num_cells + 1
    endif
    if (hmask(I+1,J  ) == 1.0 .or. hmask(I+1,J  ) == 3.0) then
      c10 = CS%bed_elev(I+1,J  ) ; num_cells = num_cells + 1
    endif
    if (hmask(I,  J+1) == 1.0 .or. hmask(I,  J+1) == 3.0) then
      c01 = CS%bed_elev(I,  J+1) ; num_cells = num_cells + 1
    endif
    if (hmask(I+1,J+1) == 1.0 .or. hmask(I+1,J+1) == 3.0) then
      c11 = CS%bed_elev(I+1,J+1) ; num_cells = num_cells + 1
    endif
    bilinear_avg = (c00 + c11) + (c10 + c01)
    if (num_cells > 0) then
      CS%bed_node(I,J) = bilinear_avg / real(num_cells)
    endif
  enddo ; enddo
  call pass_var(CS%bed_node, G%domain, position=CORNER)

  ! Step 2: Jacobi iteration to enforce cell-average constraint:
  !   (bed_node(I-1,J-1) + bed_node(I,J-1) + bed_node(I-1,J) + bed_node(I,J)) / 4 = bed_elev(i,j)
  do iter=1,max_iter
    bed_node_new(:,:) = CS%bed_node(:,:)

    do J=jsc-1,jec ; do I=isc-1,iec
      ! Node (I,J) participates in 4 cells around it. Gather each cell's residual
      ! into its own slot and reduce as (c00+c11) + (c10+c01) so the result is
      ! invariant under 90 deg rotation, which cyclically permutes the 4 cells.
      num_cells = 0
      c00 = 0.0 ; c10 = 0.0 ; c01 = 0.0 ; c11 = 0.0
      if (hmask(I,  J  ) == 1.0 .or. hmask(I,  J  ) == 3.0) then
        bilinear_avg = ((CS%bed_node(I-1,J-1) + CS%bed_node(I,  J  )) + &
                        (CS%bed_node(I,  J-1) + CS%bed_node(I-1,J  ))) * 0.25
        c00 = CS%bed_elev(I,  J  ) - bilinear_avg ; num_cells = num_cells + 1
      endif
      if (hmask(I+1,J  ) == 1.0 .or. hmask(I+1,J  ) == 3.0) then
        bilinear_avg = ((CS%bed_node(I,  J-1) + CS%bed_node(I+1,J  )) + &
                        (CS%bed_node(I+1,J-1) + CS%bed_node(I,  J  ))) * 0.25
        c10 = CS%bed_elev(I+1,J  ) - bilinear_avg ; num_cells = num_cells + 1
      endif
      if (hmask(I,  J+1) == 1.0 .or. hmask(I,  J+1) == 3.0) then
        bilinear_avg = ((CS%bed_node(I-1,J  ) + CS%bed_node(I,  J+1)) + &
                        (CS%bed_node(I,  J  ) + CS%bed_node(I-1,J+1))) * 0.25
        c01 = CS%bed_elev(I,  J+1) - bilinear_avg ; num_cells = num_cells + 1
      endif
      if (hmask(I+1,J+1) == 1.0 .or. hmask(I+1,J+1) == 3.0) then
        bilinear_avg = ((CS%bed_node(I,  J  ) + CS%bed_node(I+1,J+1)) + &
                        (CS%bed_node(I+1,J  ) + CS%bed_node(I,  J+1))) * 0.25
        c11 = CS%bed_elev(I+1,J+1) - bilinear_avg ; num_cells = num_cells + 1
      endif
      residual = (c00 + c11) + (c10 + c01)

      if (num_cells > 0) then
        bed_node_new(I,J) = CS%bed_node(I,J) + damping * residual / real(num_cells)
      endif
    enddo ; enddo

    CS%bed_node(:,:) = bed_node_new(:,:)
    call pass_var(CS%bed_node, G%domain, position=CORNER)

    ! Check convergence
    max_err = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (hmask(i,j) == 1.0 .or. hmask(i,j) == 3.0) then
        bilinear_avg = ((CS%bed_node(I-1,J-1) + CS%bed_node(I,J)) + &
                        (CS%bed_node(I,J-1)   + CS%bed_node(I-1,J))) * 0.25
        max_err = max(max_err, abs(CS%bed_elev(i,j) - bilinear_avg))
      endif
    enddo ; enddo
    call max_across_PEs(max_err)

    if (max_err < tol) exit
  enddo

  write(mesg,*) "reconstruct_bed_to_nodes max error ", max_err
  call MOM_mesg(mesg)

  if (max_err >= tol) then
    call MOM_mesg("reconstruct_bed_to_nodes: WARNING - did not converge after max_iter iterations")
  endif

end subroutine reconstruct_bed_to_nodes

!> Zero the DG(1) slope coefficients h_x, h_y at a single cell. No-op when
!! DG(1) thickness is not active. Used by external modules (e.g. melt /
!! water-flux ablation in MOM_ice_shelf) to keep the DG state self-consistent
!! when h_shelf at that cell is overwritten outside the advect step.
subroutine reset_DG_to_cellmean_at_cell(CS, i, j, h_shelf_value)
  type(ice_shelf_dyn_CS), pointer    :: CS !< Ice shelf dynamics control structure.
  integer,                intent(in) :: i  !< i index of the cell to reset.
  integer,                intent(in) :: j  !< j index of the cell to reset.
  real,                   intent(in) :: h_shelf_value !< New cell-mean thickness to broadcast to all corners [Z ~> m]

  if (.not. associated(CS)) return
  if (.not. CS%use_DG_thickness) return
  CS%h_nodal(i,j,:,:) = h_shelf_value
end subroutine reset_DG_to_cellmean_at_cell

!> Return .true. when the ice-shelf dynamics CS is associated and DG(1)
!! thickness is enabled. Lets external modules guard DG-only paths.
function is_DG_thickness_active(CS) result(active)
  type(ice_shelf_dyn_CS), pointer :: CS
  logical :: active
  active = .false.
  if (.not. associated(CS)) return
  active = CS%use_DG_thickness
end function is_DG_thickness_active

!> Reset all 4 nodal corners of every cell to zero (typical use: pre-init
!! cleanup before populating from a cell-mean field). No-op when DG(1)
!! thickness is not active.
subroutine reset_DG_to_cellmean_bulk(CS)
  type(ice_shelf_dyn_CS), pointer :: CS !< Ice shelf dynamics control structure.

  if (.not. associated(CS)) return
  if (.not. CS%use_DG_thickness) return
  CS%h_nodal(:,:,:,:) = 0.0
end subroutine reset_DG_to_cellmean_bulk

!> Accumulate a cell-mean ice-thickness source rate (positive for accumulation,
!! negative for melt) at cell (i,j) into the DG source buffer. The buffer is
!! consumed by the next ice_shelf_advect_DG1_nodal call, projected onto a
!! continuous Q1 nodal field, and applied inside the SSP-RK2 stages. No-op
!! when DG(1) thickness is not active.
subroutine accumulate_DG_source_rate(CS, i, j, rate)
  type(ice_shelf_dyn_CS), pointer    :: CS !< Ice shelf dynamics control structure.
  integer,                intent(in) :: i  !< i index of the cell.
  integer,                intent(in) :: j  !< j index of the cell.
  real,                   intent(in) :: rate !< Cell-mean thickness source rate
                                             !! to add [Z T-1 ~> m s-1].

  if (.not. associated(CS)) return
  if (.not. CS%use_DG_thickness) return
  CS%h_source_rate(i,j) = CS%h_source_rate(i,j) + rate
end subroutine accumulate_DG_source_rate

!> Project the cell-mean DG source rate CS%h_source_rate onto a continuous Q1
!! nodal source field S_node, used inside ice_shelf_advect_DG1_nodal to apply
!! basal melt + surface SMB as an RHS contribution in each SSP-RK2 stage.
!! The projection is exactly mass-conservative on hmask=1 cells: each T-cell
!! distributes its source*areaT across its 4 corners weighted by
!! CS%cell_mean_w(a,b), and each B-grid corner is the cell_mean_w-weighted
!! average of the up-to-4 hmask=1 T-cells that share it. Summing the resulting
!! corner values over a cell with the cell_mean_w weights recovers
!! areaT*S_cell exactly for cells whose corners are not at a non-hmask=1
!! boundary, and the global integral of S_node equals the global integral of
!! S_cell over hmask=1 cells. Cells with hmask=3 (Dirichlet thickness BC) are
!! excluded from the contributor pool because their h_source_rate=0 would
!! dilute the projected source at shared corners and silently lose mass from
!! the global integral. Climate-model mass deposited on hmask=3 cells is
!! accounted for separately via ISS%mass_hole in shelf_calc_flux.
!! The continuity of S_node across cell faces makes the source field smoother
!! than the per-cell S_cell, which avoids spuriously increasing DG jumps
!! where melt rates differ sharply between neighbouring cells.
subroutine project_h_source_rate_to_nodes(CS, ISS, G, S_node)
  type(ice_shelf_dyn_CS), intent(in)  :: CS  !< Ice shelf dynamics control structure.
  type(ice_shelf_state),  intent(in)  :: ISS !< Ice shelf state (hmask, h_shelf).
  type(ocean_grid_type),  intent(in)  :: G   !< The grid structure.
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(out) :: S_node !< Q1 nodal source per
                                             !! cell at the 4 corners [Z T-1 ~> m s-1].

  real, dimension(SZDIB_(G),SZDJB_(G)) :: S_corner ! Projected B-grid corner source [Z T-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: w_corner ! Total cell_mean_w summed at corner [L2]
  real :: w_contrib                                ! Per-cell-corner contribution weight [L2]
  real :: src                                      ! Cached source rate of cell (i,j) [Z T-1]
  integer :: i, j

  S_corner(:,:) = 0.0
  w_corner(:,:) = 0.0

  call pass_var(CS%h_source_rate, G%domain)
  ! Accumulate contributions from every hmask=1 T-cell to its 4 B-grid
  ! corners, weighted by CS%cell_mean_w. Only hmask=1 cells contribute, for
  ! mass conservation: cells with hmask=3 (Dirichlet thickness BC) hold
  ! h_source_rate=0 by construction, so including them in the contributor
  ! pool would add zero to num while still adding their cell_mean_w to den,
  ! diluting the projected source at shared corners and silently losing the
  ! corresponding mass from the global integral. Excluding them means each
  ! hmask=1 cell's full source is applied to its own area; mass that would
  ! be deposited on hmask=3 cells by the climate model still flows into
  ! ISS%mass_hole via the (IS_adot_int_land - adot_intt) accounting in
  ! shelf_calc_flux, unchanged from the non-DG path.
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    if (ISS%hmask(i,j) /= 1.0) cycle
    src = CS%h_source_rate(i,j)
    ! Cell-local corner (a,b) = (1,1) is the SW corner, i.e. B-node (I-1, J-1).
    w_contrib = CS%cell_mean_w(i,j,1,1)
    S_corner(I-1, J-1) = S_corner(I-1, J-1) + w_contrib * src
    w_corner(I-1, J-1) = w_corner(I-1, J-1) + w_contrib
    ! (a,b) = (2,1) = SE corner, B-node (I, J-1)
    w_contrib = CS%cell_mean_w(i,j,2,1)
    S_corner(I,   J-1) = S_corner(I,   J-1) + w_contrib * src
    w_corner(I,   J-1) = w_corner(I,   J-1) + w_contrib
    ! (a,b) = (1,2) = NW corner, B-node (I-1, J)
    w_contrib = CS%cell_mean_w(i,j,1,2)
    S_corner(I-1, J  ) = S_corner(I-1, J  ) + w_contrib * src
    w_corner(I-1, J  ) = w_corner(I-1, J  ) + w_contrib
    ! (a,b) = (2,2) = NE corner, B-node (I, J)
    w_contrib = CS%cell_mean_w(i,j,2,2)
    S_corner(I,   J  ) = S_corner(I,   J  ) + w_contrib * src
    w_corner(I,   J  ) = w_corner(I,   J  ) + w_contrib
  enddo ; enddo

  ! Normalise: at each B-grid node, S_corner becomes the cell_mean_w-weighted
  ! average of ice-covered contributing cells. Corners with no contributing
  ! ice cell get S_corner = 0 (no source there).
  do j = G%JsdB, G%JedB ; do i = G%IsdB, G%IedB
    if (w_corner(I,J) > 0.0) S_corner(I,J) = S_corner(I,J) / w_corner(I,J)
  enddo ; enddo

  ! Distribute B-grid corner values to DG cell-local corner indices for the
  ! hmask=1 cells the advect step will update. hmask=3 cells are held to
  ! their Dirichlet h_bdry_val and do not consume S_node.
  S_node(:,:,:,:) = 0.0
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle
    S_node(i,j,1,1) = S_corner(I-1, J-1)
    S_node(i,j,2,1) = S_corner(I,   J-1)
    S_node(i,j,1,2) = S_corner(I-1, J  )
    S_node(i,j,2,2) = S_corner(I,   J  )
  enddo ; enddo
end subroutine project_h_source_rate_to_nodes


! ===========================================================================
! Nodal DG(1) helpers. CS%h_nodal is the authoritative DG thickness state;
! per-cell metrics come from G%dxCv / G%dyCu / G%areaT via the Minv_xi,
! Minv_eta, cell_mean_w caches built in init_nodal_DG_metric.
! ===========================================================================

!> Read nodal-DG runtime parameters (positivity floor + driving-stress options).
subroutine read_nodal_limiter_params(param_file, mdl, CS, US)
  type(param_file_type),   intent(in)    :: param_file
  character(len=*),        intent(in)    :: mdl
  type(ice_shelf_dyn_CS),  intent(inout) :: CS
  type(unit_scale_type),   intent(in)    :: US

  call get_param(param_file, mdl, "DG1_NODAL_POSITIVITY", CS%nodal_positivity, &
                 "If true, apply the Liu-style positivity-preserving limiter to the "//&
                 "nodal DG(1) thickness corners as a safety floor against negative "//&
                 "thickness from numerical noise.", &
                 default=.true., do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG1_HIERARCHICAL_LIMITER", CS%dg_hierarchical_lim, &
                 "If true, apply a per-mode hierarchical Zhang-Shu QP-MPP slope limiter "//&
                 "with MLP-u2 vertex-based bounds (Park-Kim 2014) to the nodal DG(1) "//&
                 "thickness corners between RK stages. Decomposes the in-cell bilinear "//&
                 "polynomial into xi-slope, eta-slope, and cross-term modes and scales "//&
                 "each independently to satisfy local cell-mean bounds at the corners. "//&
                 "Mass is preserved exactly through orthogonalised mode templates against "//&
                 "cell_mean_w. References: Krivodonova 2007 JCP 226, Zhang & Shu 2010 "//&
                 "JCP 229, Park & Kim 2014 JCP 274.", &
                 default=.false., do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG1_LIMITER_ISOTROPIC", CS%dg_lim_isotropic, &
                 "If true (and DG1_HIERARCHICAL_LIMITER is also true), use the isotropic "//&
                 "single-phi Zhang-Shu MPP variant with Park-Kim MLP-u2 vertex envelopes "//&
                 "instead of the anisotropic per-mode max-product variant. The isotropic "//&
                 "form scales the full deviation (h_nodal - Hbar) at each corner by a "//&
                 "single per-cell factor; this is the textbook Zhang-Shu 2010 form and is "//&
                 "trivially rotation-symmetric and mass-conservative.", &
                 default=.false., do_not_log=(.not.CS%use_DG_thickness .or. &
                                              .not.CS%dg_hierarchical_lim))

  call get_param(param_file, mdl, "DG1_VENKAT_K", CS%dg_venkat_K, &
                 "Venkatakrishnan-style gradient-proportional slack coefficient added "//&
                 "to the surface-slope limiter envelope. Adds K*(|dS/dx|+|dS/dy|)/4 to "//&
                 "the per-cell slack, where the cell-mean S-gradient is estimated by "//&
                 "centred differences of the cell-mean surface elevation S_cell. K=1 "//&
                 "covers the corner-vs-Sbar offset of a globally linear S field exactly "//&
                 "so the limiter passes through linear gradients regardless of the "//&
                 "Park-Kim smoothness flag. K>1 allows additional curvature tolerance; "//&
                 "K=0 disables the term and reverts to Park-Kim slack only.", &
                 units="nondim", default=1.0, &
                 do_not_log=(.not.CS%use_DG_thickness .or. .not.CS%dg_hierarchical_lim))

  call get_param(param_file, mdl, "DG1_ART_VISC_C_MAX", CS%dg_art_visc_c_max, &
                 "Peak dimensionless coefficient on the DG(1) artificial-viscosity face "//&
                 "flux at fully-shocky faces. Face coefficient ramps from C_MIN (smooth) "//&
                 "to c_max (shocky) via the smoothness gate ratio r_face, which scales "//&
                 "as O(dx^2) on smooth solutions and O(1) at genuine discontinuities. "//&
                 "Smooth regions therefore receive ~0 damping, preserving optimal "//&
                 "O(dx^2) mesh convergence; only shocky faces dissipate. c_max = 0 "//&
                 "disables the viscosity. The face-node jump mode decays at 8*c_max "//&
                 "e-folds per cell traversal, so c_max = 1/16 gives one e-fold per two "//&
                 "traversals; typical 0.05-0.5. An automatic per-cell stability bound "//&
                 "(DG1_ART_VISC_KCELL) rescales all faces of any cell whose summed "//&
                 "jump-mode rates would exceed the SSP-RK2 budget.", &
                 units="nondim", default=0.0, &
                 do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG1_ART_VISC_ADVECT_COEF", CS%dg_art_visc_advect_coef, &
                 "Dimensionless multiplier on the |u_face| advective contribution to u_eff "//&
                 "in the DG(1) artificial viscosity. Per face, u_eff = advect_coef*|u_face| "//&
                 "+ strain_coef*eps_e_face*dx_perp. advect_coef = 1 (default) preserves the "//&
                 "original formulation; advect_coef = 0 drops the |u| term entirely so the "//&
                 "damping timescale is set purely by strain rate (grid-invariant by "//&
                 "construction). Useful for testing the |u| term's contribution to "//&
                 "resolution-dependent tuning.", &
                 units="nondim", default=1.0, &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))

  call get_param(param_file, mdl, "DG1_ART_VISC_STRAIN_COEF", CS%dg_art_visc_strain_coef, &
                 "Dimensionless coefficient on a velocity-independent strain-rate-scaled "//&
                 "diffusivity floor for the DG(1) artificial viscosity. Per face, "//&
                 "u_floor = strain_coef * eps_e_face * dx_perp is added to |u_face| in "//&
                 "the face flux, where eps_e_face is the 2D SSA second invariant at the "//&
                 "face midpoint. Closes the shear-margin failure mode of pure |u_face| "//&
                 "scaling (jumps excited by stretching but advectively uncoupled across "//&
                 "the face). strain_coef = 0 (default) disables the floor. Typical O(1); "//&
                 "the design rule 8*C_MAX*STRAIN_COEF = 1 sets the stagnant-region "//&
                 "damping time equal to the local strain time.", &
                 units="nondim", default=0.0, &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))

  call get_param(param_file, mdl, "DG1_ART_VISC_C_MIN", CS%dg_art_visc_c_min, &
                 "Baseline DG(1) artificial-viscosity coefficient applied at smooth "//&
                 "faces (smoothness ramp = 0). Nonzero values trade O(dx^2) smooth-"//&
                 "region convergence for damping of sub-gate jump drift.", &
                 units="nondim", default=0.0, &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))

  call get_param(param_file, mdl, "DG1_ART_VISC_WB_HARMONIC", CS%dg_art_visc_wb_harmonic, &
                 "If true, map the surface-elevation jump to the equivalent thickness "//&
                 "jump with the harmonic mean of the per-side dh/ds; if false, the "//&
                 "arithmetic mean. The two coincide on uniform-flotation faces. At "//&
                 "mixed grounded/floating faces the harmonic mean damps the surface "//&
                 "jump at the same regime-independent rate as elsewhere (m ~= 1.8 "//&
                 "rather than ~5.1 for r ~= 0.89) and makes the jump-mode stability "//&
                 "rate amplification exactly 2 at every face.", &
                 default=.true., &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))

  call get_param(param_file, mdl, "DG1_ART_VISC_GATE_SURFACE", &
                 CS%dg_art_visc_gate_surface, &
                 "If true, the DG(1) artificial-viscosity smoothness gate ratio is "//&
                 "the face surface-elevation jump relative to the local mean "//&
                 "thickness, |[s]|/H_ref, so the undamped noise floor is the same "//&
                 "size in surface elevation (the quantity that drives the spurious "//&
                 "driving-stress edge force) on grounded, floating, and mixed faces "//&
                 "alike, and DG1_ART_VISC_R_LO/R_HI read directly as the surface-"//&
                 "cliff fractions of H at which damping starts and saturates. If "//&
                 "false (legacy), the ratio is the equivalent thickness jump over "//&
                 "H_ref, which tolerates ~9x larger surface noise on grounded "//&
                 "faces than on floating ones. The R_LO/R_HI defaults depend on "//&
                 "this choice (0.01/0.1 surface mode, 0.1/1.0 legacy mode).", &
                 default=.false., &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))

  call get_param(param_file, mdl, "DG1_ART_VISC_EXCESS_JUMP", &
                 CS%dg_art_visc_excess_jump, &
                 "If true, the DG(1) artificial viscosity drives its face flux and "//&
                 "smoothness gate with only the part of the face surface-elevation "//&
                 "jump in excess of the jump supported by the adjacent cell-mean "//&
                 "surfaces. Standing discontinuity-like features that the cell "//&
                 "means also carry (e.g. shear-margin thickness contrasts) are "//&
                 "then not damped at all, while face-jump content unsupported by "//&
                 "the mean field - the spurious broken-Q1 mode - is damped at the "//&
                 "full rate. If false (legacy), the whole jump is damped, which "//&
                 "erodes persistent real contrasts at the gated rate.", &
                 default=.false., &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))
  call get_param(param_file, mdl, "DG1_ART_VISC_EXCESS_BRANCH_MAX", &
                 CS%dg_art_visc_excess_branch_max, &
                 "If true, the mean-supported allowance of DG1_ART_VISC_EXCESS_JUMP "//&
                 "is the maximum |[s]| over the admissible flotation-branch "//&
                 "assignments of each side's cell mean, where the admissible "//&
                 "branches per side are the one selected by the cell mean and the "//&
                 "one selected by the face trace (both tested at the face-QP bed). "//&
                 "This guards against under-built allowances at faces whose bed is "//&
                 "a local extremum unrepresentative of the cell interiors (e.g. a "//&
                 "grounding line sitting on a sill or ridge crest, where a floating "//&
                 "cell's mean can misread as grounded at the shallow face bed and "//&
                 "shrink the allowance, so that a legitimate margin surface step is "//&
                 "chronically damped). The allowance is unchanged wherever each "//&
                 "side's mean and trace agree on the branch, which is everywhere "//&
                 "except flotation-ambiguous faces. If false, the allowance uses "//&
                 "the mean's own branch only (legacy).", &
                 default=.true., &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0 &
                             .or. .not.CS%dg_art_visc_excess_jump))

  ! NOTE: DG1_ART_VISC_GATE_SURFACE must be read before R_LO/R_HI: their defaults
  ! depend on the gate mode.
  call get_param(param_file, mdl, "DG1_ART_VISC_R_LO", CS%dg_art_visc_r_lo, &
                 "Lower threshold of the DG(1) artificial-viscosity smoothness "//&
                 "gate: faces with gate ratio r_face below this receive the "//&
                 "baseline coefficient only. In surface-gate mode "//&
                 "(DG1_ART_VISC_GATE_SURFACE) r_face = |[s]|/H_ref, so R_LO is "//&
                 "the surface-cliff height, as a fraction of the local mean "//&
                 "thickness, at which damping starts; in legacy mode r_face = "//&
                 "|Delta h_eq|/H_ref (a relative thickness jump). Under "//&
                 "DG1_ART_VISC_EXCESS_JUMP the jump in either ratio is the mean-"//&
                 "unsupported excess. Either ratio is "//&
                 "O(dx^2) on smooth solutions, so a fixed threshold classifies "//&
                 "smoothness resolution-invariantly; jumps below the threshold "//&
                 "persist undamped.", &
                 units="nondim", default=merge(0.01, 0.1, CS%dg_art_visc_gate_surface), &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))
  call get_param(param_file, mdl, "DG1_ART_VISC_R_HI", CS%dg_art_visc_r_hi, &
                 "Saturation threshold of the DG(1) artificial-viscosity smoothness "//&
                 "gate: faces with gate ratio at or above this receive the full "//&
                 "C_MAX. See DG1_ART_VISC_R_LO for the gate-mode-dependent "//&
                 "meaning of the ratio.", &
                 units="nondim", default=merge(0.1, 1.0, CS%dg_art_visc_gate_surface), &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))
  call get_param(param_file, mdl, "DG1_ART_VISC_KCELL", CS%dg_art_visc_kcell, &
                 "Per-cell stability budget for the DG(1) artificial viscosity: the "//&
                 "sum over a cell's faces of the jump-mode decay rates "//&
                 "(4*amp*c*u_eff/dx_perp) times dt is limited to this value by "//&
                 "rescaling the cell's face coefficients. SSP-RK2 requires < 2; "//&
                 "the default 1.0 leaves a 2x margin.", &
                 units="nondim", default=1.0, &
                 do_not_log=(.not.CS%use_DG_thickness .or. CS%dg_art_visc_c_max == 0.0))
  if ((CS%use_DG_thickness) .and. (CS%dg_art_visc_c_max > 0.0)) then
    if (CS%dg_art_visc_r_hi <= CS%dg_art_visc_r_lo) call MOM_error(FATAL, &
        "MOM_ice_shelf_dynamics, initialize_ice_shelf_dyn: DG1_ART_VISC_R_HI must "//&
        "exceed DG1_ART_VISC_R_LO.")
  endif

  ! Stagnant-jump diagnostic thresholds. Internal constants, not user knobs.
  CS%dg_slow_idle_u_tiny   = 1.0e-9  * US%m_s_to_L_T
  CS%dg_slow_idle_eps_tiny = 1.0e-12 * US%s_to_T
  CS%dg_slow_idle_s_tol    = 1.0     * US%m_to_Z

  call get_param(param_file, mdl, "DG_DRIVING_STRESS_IBP", CS%dg_driving_stress_IBP, &
                 "If true, evaluate the DG(1) driving stress with the integration-by-parts "//&
                 "weak form. Interior faces use a central numerical flux "//&
                 "P_star = 0.5*(P_loc + P_ngh); any single-valued P_star gives the same SSA "//&
                 "node-assembly total, so the prior Rusanov/IIPG penalty machinery is gone. "//&
                 "If false (default), evaluate -rho*g*h*grad(s) directly at 2x2 Gauss points "//&
                 "using the Q1 nodal basis, with an optional scale-aware face flux gated by "//&
                 "DG_FACE_FLUX_K_THRESH.", &
                 default=.false., do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG_FACE_FLUX_K_THRESH", CS%dg_face_flux_K_thresh, &
                 "Controls the optional face flux at interior hmask=1 / hmask=1 faces in "//&
                 "the strong-form DG(1) driving stress. K_THRESH < 0 (the default) disables "//&
                 "the face flux entirely (pure strong form). K_THRESH = 0 forces s = 1 "//&
                 "(full central-IBP face flux: equivalent at the SSA node assembly to the "//&
                 "IBP routine with any single-valued P_star). K_THRESH > 0 enables a "//&
                 "Venkatakrishnan-style smooth blend s = r^2 / (r^2 + K^2) with "//&
                 "r = |[h]| / h_avg, ramping the per-cell face contribution from zero at "//&
                 "smooth faces to the central-IBP face flux at large jumps. Typical "//&
                 "engaging values are ~0.05-0.2. Values >> 1 are effectively disabled "//&
                 "because r << 1 in typical near-C0 nodal Q1 flows. Has no effect when "//&
                 "USE_DG_THICKNESS is false or DG_DRIVING_STRESS_IBP is true.", &
                 units="nondim", default=-1.0, do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG_GL_GATE_CONTINUOUS", CS%dg_gl_gate_continuous, &
                 "If true, evaluate the basal-friction flotation gates, ground_frac, "//&
                 "and the Coulomb effective pressure on a continuous node-averaged "//&
                 "thickness field instead of each cell's own DG nodal thickness. This "//&
                 "prevents the flotation crossing from hiding inside an inter-cell "//&
                 "thickness jump, which otherwise gives the basal friction a dead band "//&
                 "(no response to thickness changes while the crossing traverses the "//&
                 "jump) and lets the grounding line lock at cell faces. Force "//&
                 "magnitudes always use the true DG thickness, and the driving-stress "//&
                 "branch selectors keep per-side own-thickness tests unless "//&
                 "DG_GL_GATE_DRIVING_STRESS is also true. Has no effect on the IBP "//&
                 "driving-stress path (DG_DRIVING_STRESS_IBP = True) beyond its "//&
                 "influence on ground_frac.", &
                 default=.false., do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG_GL_GATE_DRIVING_STRESS", CS%dg_gl_gate_driving_stress, &
                 "If true (with DG_GL_GATE_CONTINUOUS), the strong-form driving-stress "//&
                 "flotation branches (volume, subgrid, and face Dirac term) also use the "//&
                 "continuous gate field. If false (default), the driving stress keeps "//&
                 "per-side own-thickness branching, equivalent to the continuous "//&
                 "s = max(h - bed, (1-rho_i/rho_w)*h), so the assembled force is "//&
                 "continuous in the model state. Gating the face Dirac term with the "//&
                 "continuous field makes its [s] switch discontinuously between [h] and "//&
                 "(1-rho_i/rho_w)*[h] when the gate sign flips at a face, a force "//&
                 "discontinuity of order half a cell's driving stress that pins the "//&
                 "steady grounding line at cell faces (where the node-average deficit "//&
                 "is zero). Only the basal-friction gates, ground_frac, and Coulomb "//&
                 "effective pressure should normally use the continuous gate.", &
                 default=.false., do_not_log=.not.(CS%use_DG_thickness .and. CS%dg_gl_gate_continuous))

  call get_param(param_file, mdl, "DG_GL_GATE_DEFICIT_SCALE", CS%dg_gl_gate_deficit_scale, &
                 "Regularization scale s for inverse-flotation-deficit weighting of the "//&
                 "DG_GL_GATE_CONTINUOUS node average. <= 0 gives the plain arithmetic "//&
                 "corner mean. > 0 weights each corner of the gate field by "//&
                 "1/(|(rho_i/rho_w)*h - bed| + s), so the side nearer flotation "//&
                 "dominates and a large one-sided thickness jump at the grounding-line "//&
                 "face cannot drag the gate's flotation crossing far into the lighter "//&
                 "cell. Small s approaches pinning straddling faces exactly at "//&
                 "flotation, which makes the gate insensitive to thickness changes "//&
                 "while the face straddles (a dead band) and should be avoided; large "//&
                 "s approaches the arithmetic mean. A few meters is a reasonable "//&
                 "starting value. Meaningful only when DG_GL_GATE_CONTINUOUS is true.", &
                 units="m", default=0.0, scale=US%m_to_Z, &
                 do_not_log=.not.(CS%use_DG_thickness .and. CS%dg_gl_gate_continuous))

end subroutine read_nodal_limiter_params

!> Initialise the per-cell metric tables (Minv_xi, Minv_eta, cell_mean_w) from
!! the grid. Per-cell face lengths come from G%dxCv (south/north) and G%dyCu
!! (west/east); the bilinear face-length interpolation is
!!   a(eta) = dxS*(1-eta) + dxN*eta   on eta in [0,1]
!!   d(xi)  = dyW*(1-xi)  + dyE*xi    on xi  in [0,1]
!! with face-length aliases dxS = G%dxCv(i,J-1), dxN = G%dxCv(i,J),
!! dyW = G%dyCu(I-1,j), dyE = G%dyCu(I,j). Mass-matrix factors are analytic for
!! the locally-orthogonal lat/lon grid: M = M_xi (x) M_eta with
!! M_xi_{a,a'}  = int N_a(xi)*N_a'(xi)*d(xi) dxi
!! M_eta_{b,b'} = int N_b(eta)*N_b'(eta)*a(eta) deta,
!! where N_1 = 1-x, N_2 = x. The 2x2 inverses are analytic.
subroutine init_nodal_DG_metric(CS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ocean_grid_type),  intent(in)    :: G

  real :: dxS, dxN, dyW, dyE    ! face lengths [L ~> m]
  real :: M11, M12, M22, det    ! mass-matrix entries and determinant
  integer :: i, j, isd, ied, jsd, jed

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j = jsd, jed ; do i = isd, ied
    ! Face lengths with one-sided fallback at non-reentrant domain edges. For
    ! reentrant domains the wrap halo carries valid dxCv/dyCu, so use the
    ! two-sided metric even at the global west/south edge cell.
    if ((J-1 >= G%JsdB) .and. (CS%reentrant_y .or. (j + G%jdg_offset > G%jsg))) then
      dxS = G%dxCv(i,J-1) ; dxN = G%dxCv(i,J)
    else
      dxS = G%dxCv(i,J)   ; dxN = G%dxCv(i,J)
    endif
    if ((I-1 >= G%IsdB) .and. (CS%reentrant_x .or. (i + G%idg_offset > G%isg))) then
      dyW = G%dyCu(I-1,j) ; dyE = G%dyCu(I,j)
    else
      dyW = G%dyCu(I,j)   ; dyE = G%dyCu(I,j)
    endif

    ! M_xi: int_0^1 N_a*N_a'*d(xi) dxi where d(xi) = dyW*(1-xi) + dyE*xi.
    ! Closed-form: M11 = dyW/4 + dyE/12, M22 = dyW/12 + dyE/4, M12 = dyW/12 + dyE/12.
    M11 = dyW/4.0 + dyE/12.0
    M22 = dyW/12.0 + dyE/4.0
    M12 = (dyW + dyE)/12.0
    det = M11*M22 - M12*M12
    if (det > 0.0) then
      CS%Minv_xi(i,j,1,1) =  M22 / det
      CS%Minv_xi(i,j,2,2) =  M11 / det
      CS%Minv_xi(i,j,1,2) = -M12 / det
      CS%Minv_xi(i,j,2,1) = -M12 / det
    endif

    ! M_eta: int_0^1 N_b*N_b'*a(eta) deta where a(eta) = dxS*(1-eta) + dxN*eta.
    M11 = dxS/4.0 + dxN/12.0
    M22 = dxS/12.0 + dxN/4.0
    M12 = (dxS + dxN)/12.0
    det = M11*M22 - M12*M12
    if (det > 0.0) then
      CS%Minv_eta(i,j,1,1) =  M22 / det
      CS%Minv_eta(i,j,2,2) =  M11 / det
      CS%Minv_eta(i,j,1,2) = -M12 / det
      CS%Minv_eta(i,j,2,1) = -M12 / det
    endif

    ! Per-corner integration weight w(a,b) = int_0^1 int_0^1 N(a,b)*a(eta)*d(xi) dxi deta.
    ! N(1,1) = (1-xi)(1-eta), etc. Closed-form by separation:
    !   w(a,b) = (int N_a(xi)*d(xi) dxi) * (int N_b(eta)*a(eta) deta)
    ! int N_1(xi)*d(xi) dxi = dyW/3 + dyE/6
    ! int N_2(xi)*d(xi) dxi = dyW/6 + dyE/3
    ! int N_1(eta)*a(eta) deta = dxS/3 + dxN/6
    ! int N_2(eta)*a(eta) deta = dxS/6 + dxN/3
    CS%cell_mean_w(i,j,1,1) = (dyW/3.0 + dyE/6.0) * (dxS/3.0 + dxN/6.0)
    CS%cell_mean_w(i,j,2,1) = (dyW/6.0 + dyE/3.0) * (dxS/3.0 + dxN/6.0)
    CS%cell_mean_w(i,j,1,2) = (dyW/3.0 + dyE/6.0) * (dxS/6.0 + dxN/3.0)
    CS%cell_mean_w(i,j,2,2) = (dyW/6.0 + dyE/3.0) * (dxS/6.0 + dxN/3.0)
  enddo ; enddo

  ! Cache the orthogonalisation offsets for the hierarchical-limiter mode
  ! templates. The standard polynomial modes (xi-0.5, eta-0.5, (xi-0.5)(eta-0.5))
  ! are not exactly orthogonal to the constant under non-uniform cell_mean_w
  ! (curvilinear cells); subtracting their weighted-mean offset gives modes
  ! that are exactly orthogonal to the constant, so per-mode scaling preserves
  ! Hbar exactly. Offsets are tiny on near-Cartesian grids and zero on
  ! perfectly uniform ones. Stored once at init since cell_mean_w is fixed.
  if (allocated(CS%mu_lim_xi)) then
    do j = jsd, jed ; do i = isd, ied
      call compute_mu_lim_offsets(CS%cell_mean_w(i,j,:,:), &
                                  CS%mu_lim_xi(i,j), CS%mu_lim_eta(i,j), CS%mu_lim_cross(i,j))
    enddo ; enddo
  endif
end subroutine init_nodal_DG_metric

!> Helper: orthogonalisation offsets so that the modal templates
!! B_orth, C_orth, D_orth are exactly cell_mean_w-orthogonal to the constant
!! mode (i.e., have zero weighted mean over the cell corners).
pure subroutine compute_mu_lim_offsets(w_cell, mu_B, mu_C, mu_D)
  real, dimension(2,2), intent(in)  :: w_cell  !< Per-cell cell_mean_w pre-normalised weights
  real,                 intent(out) :: mu_B, mu_C, mu_D
  real :: area
  area = ((w_cell(1,1) + w_cell(2,2)) + (w_cell(1,2) + w_cell(2,1)))
  if (area > 0.0) then
    mu_B = (((-0.5)*w_cell(1,1) + (0.5)*w_cell(2,1)) + &
            ((-0.5)*w_cell(1,2) + (0.5)*w_cell(2,2))) / area
    mu_C = (((-0.5)*w_cell(1,1) + (-0.5)*w_cell(2,1)) + &
            (( 0.5)*w_cell(1,2) + ( 0.5)*w_cell(2,2))) / area
    mu_D = ((( 0.25)*w_cell(1,1) + (-0.25)*w_cell(2,1)) + &
            ((-0.25)*w_cell(1,2) + ( 0.25)*w_cell(2,2))) / area
  else
    mu_B = 0.0 ; mu_C = 0.0 ; mu_D = 0.0
  endif
end subroutine compute_mu_lim_offsets

!> For symmetric BGRID + reentrant domains, the wrap-mate B-nodes on either
!! side of the periodic boundary are the same physical location, but each
!! cell adjacent to the boundary stores its own DG corner there. At IC we
!! want the two wrap-shared corners to hold the same value (a continuous
!! initial field). FMS pass_var with position=CORNER does not always
!! reconcile owned wrap-edge B-nodes for symmetric grids, and the IC node
!! file may have small asymmetries at the wrap. Force consistency by
!! having the PE owning the global west (south) edge adopt the wrap-mate's
!! value from its halo. Halos are then re-passed.
subroutine enforce_wrap_corner_consistency(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ice_shelf_state),  intent(in)    :: ISS
  type(ocean_grid_type),  intent(inout) :: G

  integer :: i, j

  if (.not. G%symmetric) return
  if (.not. (CS%reentrant_x .or. CS%reentrant_y)) return

  if (CS%reentrant_x .and. (G%isc + G%idg_offset == G%isg)) then
    i = G%isc
    do j = G%jsc, G%jec
      if (ISS%hmask(i, j) == 1.0 .or. ISS%hmask(i, j) == 3.0) then
        CS%h_nodal(i, j, 1, 1) = CS%h_nodal(i-1, j, 2, 1)
        CS%h_nodal(i, j, 1, 2) = CS%h_nodal(i-1, j, 2, 2)
      endif
    enddo
  endif

  if (CS%reentrant_y .and. (G%jsc + G%jdg_offset == G%jsg)) then
    j = G%jsc
    do i = G%isc, G%iec
      if (ISS%hmask(i, j) == 1.0 .or. ISS%hmask(i, j) == 3.0) then
        CS%h_nodal(i, j, 1, 1) = CS%h_nodal(i, j-1, 1, 2)
        CS%h_nodal(i, j, 2, 1) = CS%h_nodal(i, j-1, 2, 2)
      endif
    enddo
  endif

  ! pass_corner_field is collective; every PE must call regardless of
  ! whether it touched its data.
  call pass_corner_field(CS%h_nodal, G)
end subroutine enforce_wrap_corner_consistency

!> Halo-exchange the per-cell 4-corner nodal field. Each corner slot is a
!! cell-centered scalar (A-grid).
subroutine pass_corner_field(h_nodal, G)
  type(ocean_grid_type),  intent(in) :: G
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(inout) :: h_nodal
  real, dimension(SZDI_(G),SZDJ_(G),4) :: tmp

  tmp(:,:,1) = h_nodal(:,:,1,1)
  tmp(:,:,2) = h_nodal(:,:,2,1)
  tmp(:,:,3) = h_nodal(:,:,1,2)
  tmp(:,:,4) = h_nodal(:,:,2,2)
  call pass_var(tmp, G%domain)
  h_nodal(:,:,1,1) = tmp(:,:,1)
  h_nodal(:,:,2,1) = tmp(:,:,2)
  h_nodal(:,:,1,2) = tmp(:,:,3)
  h_nodal(:,:,2,2) = tmp(:,:,4)
end subroutine pass_corner_field

!> Compute the area-weighted cell mean of the 4 corner values via cell_mean_w.
pure real function nodal_cell_mean(h_cell, w_cell) result(Hbar)
  real, dimension(2,2), intent(in) :: h_cell, w_cell
  real :: area
  area = ((w_cell(1,1) + w_cell(2,2)) + (w_cell(1,2) + w_cell(2,1)))
  if (area > 0.0) then
    Hbar = ( (w_cell(1,1)*h_cell(1,1) + w_cell(2,2)*h_cell(2,2)) + &
             (w_cell(1,2)*h_cell(1,2) + w_cell(2,1)*h_cell(2,1)) ) / area
  else
    Hbar = 0.0
  endif
end function nodal_cell_mean

!> Publish ISS%h_shelf from CS%h_nodal as the area-weighted mean.
subroutine recompute_h_shelf_from_nodal(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(in) :: CS
  type(ice_shelf_state),  intent(in) :: ISS  ! h_shelf is a pointer; writing through it is OK
  type(ocean_grid_type),  intent(in) :: G

  integer :: i, j
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    if (ISS%hmask(i,j) == 1.0 .or. ISS%hmask(i,j) == 3.0) then
      ISS%h_shelf(i,j) = nodal_cell_mean(CS%h_nodal(i,j,:,:), CS%cell_mean_w(i,j,:,:))
    endif
  enddo ; enddo
end subroutine recompute_h_shelf_from_nodal

!> Cold-start initialise h_nodal from the cell-mean h_shelf field. Each
!! corner is the area-weighted average of the up-to-four surrounding cells'
!! h_shelf values, using G%areaT as the metric weight. On uniform grids this
!! collapses bit-identically to the simple arithmetic mean. One-sided
!! fallback at hmask boundaries. Produces a continuous (no-jump) Q1 field at
!! init; the limiter introduces admissible jumps later if the field warrants.
subroutine initialize_h_nodal_from_cellmean(h_shelf, h_nodal, hmask, G)
  type(ocean_grid_type), intent(inout) :: G
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in) :: h_shelf
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(inout) :: h_nodal
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in) :: hmask

  integer :: i, j, a, b, ic, jc, dx, dy, di_off, dj_off
  real :: sum_area_h ! Area-weighted sum of h_shelf values [Z L2 ~> m3]
  real :: sum_area   ! Sum of areas of contributing cells [L2 ~> m2]
  real :: area_ic    ! Area of one contributing cell [L2 ~> m2]

  ! Corner (a,b) of cell (i,j) is shared with up to four T-cells whose
  ! offsets are (dx*di_off, dy*dj_off) for dx,dy in {0,1} and
  ! di_off = 2*(a-1)-1, dj_off = 2*(b-1)-1 (i.e. -1 for a=1, +1 for a=2).
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (hmask(i,j) /= 1.0 .and. hmask(i,j) /= 3.0) cycle
    do b = 1, 2 ; do a = 1, 2
      di_off = 2*(a-1) - 1
      dj_off = 2*(b-1) - 1
      sum_area_h = 0.0 ; sum_area = 0.0
      do dy = 0, 1 ; do dx = 0, 1
        ic = i + dx*di_off
        jc = j + dy*dj_off
        if (ic < G%isd .or. ic > G%ied) cycle
        if (jc < G%jsd .or. jc > G%jed) cycle
        if (hmask(ic,jc) == 1.0 .or. hmask(ic,jc) == 3.0) then
          area_ic = G%areaT(ic, jc)
          sum_area_h = sum_area_h + area_ic * h_shelf(ic, jc)
          sum_area   = sum_area   + area_ic
        endif
      enddo ; enddo
      if (sum_area > 0.0) then
        h_nodal(i,j,a,b) = sum_area_h / sum_area
      else
        h_nodal(i,j,a,b) = h_shelf(i,j)
      endif
    enddo ; enddo
  enddo ; enddo
end subroutine initialize_h_nodal_from_cellmean


!> Liu-style positivity-preserving limiter: scale each cell's corner
!! deviations from the mean by a single factor in [0,1] so the minimum
!! corner is at least CS%min_h_shelf. Preserves cell mean exactly.
subroutine nodal_positivity_limit(CS, G, ISS)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ocean_grid_type),  intent(inout) :: G
  type(ice_shelf_state),  intent(in)    :: ISS

  integer :: i, j, a, b
  real :: Hbar, hmin, phi_pos, denom

  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle
    Hbar = nodal_cell_mean(CS%h_nodal(i,j,:,:), CS%cell_mean_w(i,j,:,:))
    hmin = min(min(CS%h_nodal(i,j,1,1), CS%h_nodal(i,j,2,2)), &
               min(CS%h_nodal(i,j,1,2), CS%h_nodal(i,j,2,1)))
    if (hmin >= CS%min_h_shelf) cycle
    denom = Hbar - hmin
    if (denom > 1.0e-30) then
      phi_pos = max(0.0, min(1.0, (Hbar - CS%min_h_shelf)/denom))
    else
      phi_pos = 0.0
    endif
    do b = 1, 2 ; do a = 1, 2
      CS%h_nodal(i,j,a,b) = Hbar + phi_pos*(CS%h_nodal(i,j,a,b) - Hbar)
    enddo ; enddo
  enddo ; enddo

  call pass_corner_field(CS%h_nodal, G)
end subroutine nodal_positivity_limit

!> Per-mode hierarchical Zhang-Shu QP-MPP slope limiter for DG(1) thickness
!! with MLP-u2 vertex-based bounds. Decomposes each cell's in-cell bilinear
!! polynomial into (xi-slope, eta-slope, cross-term) modes and scales each
!! independently so the limited corner values stay within the per-corner
!! local cell-mean envelope drawn from the cells touching each B-node
!! (Park-Kim 2014 vertex-based MLP family). Limiting hierarchy: cross first,
!! then xi-slope, then eta-slope, each pass using the already-limited
!! contributions from previous passes as part of the residual budget.
!!
!! Mass is preserved exactly: the mode templates are pre-orthogonalised
!! against the cell_mean_w-weighted constant mode (cached as mu_lim_*),
!! so each mode's contribution has zero weighted mean by construction.
!!
!! References: Krivodonova 2007 JCP 226 (hierarchical mode-by-mode limiting),
!! Zhang & Shu 2010 JCP 229 (QP-MPP theorem), Park & Kim 2014 JCP 274
!! (MLP-u2 vertex-based bound construction).
subroutine nodal_hierarchical_limit(CS, G, ISS)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ocean_grid_type),  intent(in)    :: G
  type(ice_shelf_state),  intent(in)    :: ISS

  real, dimension(SZDIB_(G),SZDJB_(G)) :: Hmax_B, Hmin_B
  real, dimension(SZDIB_(G),SZDJB_(G)) :: count_B
  real, dimension(SZDI_(G),SZDJ_(G))   :: H_cell    ! Cell-mean field for Park-Kim stencil [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G))   :: pk_factor ! Park-Kim smooth-extrema factor [nondim, 0..1]
  real, dimension(SZDI_(G),SZDJ_(G))   :: pk_d2x_cell, pk_d2y_cell  ! Cell-mean second diffs [Z ~> m]
  real :: pk_slack
  real :: phi_iso, iso_dev                  ! Isotropic-variant scratch
  real, parameter :: H_LARGE = 1.0e30
  real, parameter :: TINY_DEV = 1.0e-30    ! Threshold for skipping near-zero deviations
  ! Park-Kim slack coefficients. PK_SLACK_D2 is the natural-overshoot
  ! coefficient for a smooth Q1 quadratic peak: corner - bound ~ |d2|/12.
  ! PK_SLACK_ENV is a small extra envelope-width margin to absorb
  ! sub-quadratic curvature. Both apply only at smooth-flagged cells.
  real, parameter :: PK_SLACK_D2  = 1.0/12.0  ! ~0.083
  real, parameter :: PK_SLACK_ENV = 0.05      ! 5% of envelope width
  ! Bound-tolerance constants: relax the MLP-u2 vertex envelope by
  ! max(BOUND_TOL_ABS, BOUND_TOL_REL * envelope_width) before applying
  ! the per-mode MPP scaling. Prevents spurious limiting at smooth local
  ! extrema of the cell-mean field where the home cell happens to be the
  ! envelope extremum at one corner. Both constants are well below any
  ! physically meaningful overshoot for ice-shelf thickness (~1 m).
  real, parameter :: BOUND_TOL_ABS = 1.0e-4  ! [Z ~> m] absolute floor
  real, parameter :: BOUND_TOL_REL = 1.0e-6  ! [nondim] relative
  ! Tolerances for the max-product iterative projection.
  real, parameter :: LP_TOL = 1.0e-10      ! [nondim] feasibility tolerance
  integer, parameter :: MP_MAX_ITER = 20   ! Iteration cap for max-product convergence
  real :: cell_mean_val
  real :: Hbar, Hbar_new
  real :: h11, h21, h12, h22
  real :: bb, cc, dd                         ! modal coefficients (uniform-w formulas)
  real :: muB, muC, muD                      ! orthogonalisation offsets for this cell
  real :: B_orth(2,2), C_orth(2,2), D_orth(2,2)
  real :: dev_b(2,2), dev_c(2,2), dev_d(2,2) ! per-corner mode contributions
  real :: bound_max(2,2), bound_min(2,2)
  real :: bound_tol                          ! Per-corner relaxation of the envelope
  real :: phi_b, phi_c, phi_d
  ! Max-product scratch variables (constraints are A(:,k) . phi <= B(k) for k = 1..8).
  real :: lp_A(3,8), lp_B(8)
  real :: lp_inv_det, lp_x(3), lp_best_obj, lp_obj
  integer :: lp_i, lp_k, lp_q
  integer :: mp_n_pos                       ! Number of positive-coef modes at the active constraint
  real :: budget_high, budget_low, ratio
  real :: h_new(2,2)
  integer :: i, j, a, b, I_node, J_node

  if (.not. associated(CS%h_nodal)) return

  ! Build vertex-based cell-mean envelope at each B-node (Park-Kim MLP-u2
  ! style). Includes hmask=1 and hmask=3 (Dirichlet) cells; skips hmask=0/2.
  call pass_corner_field(CS%h_nodal, G)
  Hmax_B(:,:)  = -H_LARGE
  Hmin_B(:,:)  =  H_LARGE
  count_B(:,:) = 0.0
  ! Also build a 2D cell-mean field H_cell for the Park-Kim smooth-extrema
  ! indicator below. Sentinel value H_LARGE marks unavailable cells (hmask=0
  ! or hmask=2). Both the envelope and H_cell are built in the same pass.
  H_cell(:,:) = H_LARGE
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    if (ISS%hmask(i,j) == 1.0) then
      cell_mean_val = nodal_cell_mean(CS%h_nodal(i,j,:,:), CS%cell_mean_w(i,j,:,:))
    elseif (ISS%hmask(i,j) == 3.0) then
      cell_mean_val = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
    else
      cycle
    endif
    H_cell(i, j) = cell_mean_val
    if (i-1 >= G%IsdB .and. j-1 >= G%JsdB) then
      Hmax_B(i-1, j-1) = max(Hmax_B(i-1, j-1), cell_mean_val)
      Hmin_B(i-1, j-1) = min(Hmin_B(i-1, j-1), cell_mean_val)
      count_B(i-1, j-1) = count_B(i-1, j-1) + 1.0
    endif
    if (i   <= G%IedB .and. j-1 >= G%JsdB) then
      Hmax_B(i,   j-1) = max(Hmax_B(i,   j-1), cell_mean_val)
      Hmin_B(i,   j-1) = min(Hmin_B(i,   j-1), cell_mean_val)
      count_B(i,   j-1) = count_B(i,   j-1) + 1.0
    endif
    if (i-1 >= G%IsdB .and. j   <= G%JedB) then
      Hmax_B(i-1, j  ) = max(Hmax_B(i-1, j  ), cell_mean_val)
      Hmin_B(i-1, j  ) = min(Hmin_B(i-1, j  ), cell_mean_val)
      count_B(i-1, j  ) = count_B(i-1, j  ) + 1.0
    endif
    if (i   <= G%IedB .and. j   <= G%JedB) then
      Hmax_B(i,   j  ) = max(Hmax_B(i,   j  ), cell_mean_val)
      Hmin_B(i,   j  ) = min(Hmin_B(i,   j  ), cell_mean_val)
      count_B(i,   j  ) = count_B(i,   j  ) + 1.0
    endif
  enddo ; enddo
  call pass_var(Hmax_B,  G%domain, position=CORNER)
  call pass_var(Hmin_B,  G%domain, position=CORNER)
  call pass_var(count_B, G%domain, position=CORNER)

  ! Park-Kim MLP-u2 smooth-extrema indicator (Park & Kim 2014, JCP 274).
  ! Uses second-difference sign consistency across a wider stencil: a cell
  ! is at a smooth extremum if its second-difference d2 has the same sign
  ! as its neighbors' d2 in each direction. Smooth peaks have consistently
  ! negative d2 over a 3-cell neighborhood; oscillations have alternating
  ! signs. The sign-consistency check requires d2 at neighbouring cells,
  ! hence Hbar in a 2-cell-deep halo. A pass_var on H_cell extends the
  ! field across PE boundaries to the required depth.
  call pass_var(H_cell, G%domain)
  ! Step 1: compute d2x, d2y at every cell where the immediate-neighbour
  ! cell means are available (1-deep halo from the owned cells).
  pk_d2x_cell(:,:) = 0.0
  pk_d2y_cell(:,:) = 0.0
  do j = G%jsd+1, G%jed-1 ; do i = G%isd+1, G%ied-1
    if (H_cell(i, j) >= H_LARGE - 1.0) cycle
    if (H_cell(i-1, j) < H_LARGE - 1.0 .and. H_cell(i+1, j) < H_LARGE - 1.0) then
      pk_d2x_cell(i, j) = (H_cell(i-1, j) - 2.0 * H_cell(i, j)) + H_cell(i+1, j)
    endif
    if (H_cell(i, j-1) < H_LARGE - 1.0 .and. H_cell(i, j+1) < H_LARGE - 1.0) then
      pk_d2y_cell(i, j) = (H_cell(i, j-1) - 2.0 * H_cell(i, j)) + H_cell(i, j+1)
    endif
  enddo ; enddo

  ! Step 2: per owned cell, check sign consistency of d2 in x and y over
  ! the 3-cell stencil centred at K. If consistent in both directions, the
  ! cell is at a smooth extremum and pk_factor = 1; otherwise 0.
  pk_factor(:,:) = 0.0
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle
    if (i-1 < G%isd+1 .or. i+1 > G%ied-1) cycle
    if (j-1 < G%jsd+1 .or. j+1 > G%jed-1) cycle
    ! d2 same sign in x at i-1, i, i+1?
    if (((pk_d2x_cell(i-1, j) >= 0.0) .and. (pk_d2x_cell(i, j) >= 0.0) .and. &
         (pk_d2x_cell(i+1, j) >= 0.0)) .or. &
        ((pk_d2x_cell(i-1, j) <= 0.0) .and. (pk_d2x_cell(i, j) <= 0.0) .and. &
         (pk_d2x_cell(i+1, j) <= 0.0))) then
      ! d2 same sign in y at j-1, j, j+1?
      if (((pk_d2y_cell(i, j-1) >= 0.0) .and. (pk_d2y_cell(i, j) >= 0.0) .and. &
           (pk_d2y_cell(i, j+1) >= 0.0)) .or. &
          ((pk_d2y_cell(i, j-1) <= 0.0) .and. (pk_d2y_cell(i, j) <= 0.0) .and. &
           (pk_d2y_cell(i, j+1) <= 0.0))) then
        pk_factor(i, j) = 1.0
      endif
    endif
  enddo ; enddo

  ! Reset diagnostic buffers (limiter inactive at non-hmask=1 cells -> phi=1, drift=0).
  if (associated(CS%dg_lim_phi_xi))     CS%dg_lim_phi_xi(:,:)     = 1.0
  if (associated(CS%dg_lim_phi_eta))    CS%dg_lim_phi_eta(:,:)    = 1.0
  if (associated(CS%dg_lim_phi_cross))  CS%dg_lim_phi_cross(:,:)  = 1.0
  if (associated(CS%dg_lim_mass_drift)) CS%dg_lim_mass_drift(:,:) = 0.0
  if (associated(CS%dg_lim_phi))        CS%dg_lim_phi(:,:)        = 1.0
  if (associated(CS%dg_lim_pk_factor))  CS%dg_lim_pk_factor(:,:)  = pk_factor(:,:)

  ! Per-cell hierarchical limiting.
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle

    ! Cell mean (exact, using non-uniform cell_mean_w).
    Hbar = nodal_cell_mean(CS%h_nodal(i,j,:,:), CS%cell_mean_w(i,j,:,:))

    h11 = CS%h_nodal(i,j,1,1) ; h21 = CS%h_nodal(i,j,2,1)
    h12 = CS%h_nodal(i,j,1,2) ; h22 = CS%h_nodal(i,j,2,2)

    ! Standard polynomial mode coefficients (uniform-weight formulas). Under
    ! non-uniform cell_mean_w these are not the exact L2(w) projections but
    ! the *reconstruction* below uses orthogonalised templates that preserve
    ! Hbar exactly, so mass conservation is unaffected.
    bb = 0.5 * ((h21 + h22) - (h11 + h12))            ! east-mean - west-mean
    cc = 0.5 * ((h12 + h22) - (h11 + h21))            ! north-mean - south-mean
    dd = (h11 + h22) - (h12 + h21)                    ! diagonal vs anti-diagonal

    ! Orthogonalised mode templates per corner: standard polynomial values
    ! minus their cell_mean_w-weighted means (cached).
    muB = CS%mu_lim_xi(i,j) ; muC = CS%mu_lim_eta(i,j) ; muD = CS%mu_lim_cross(i,j)
    B_orth(1,1) = -0.5  - muB ; B_orth(2,1) =  0.5  - muB
    B_orth(1,2) = -0.5  - muB ; B_orth(2,2) =  0.5  - muB
    C_orth(1,1) = -0.5  - muC ; C_orth(2,1) = -0.5  - muC
    C_orth(1,2) =  0.5  - muC ; C_orth(2,2) =  0.5  - muC
    D_orth(1,1) =  0.25 - muD ; D_orth(2,1) = -0.25 - muD
    D_orth(1,2) = -0.25 - muD ; D_orth(2,2) =  0.25 - muD

    do b = 1, 2 ; do a = 1, 2
      dev_b(a,b) = bb * B_orth(a,b)
      dev_c(a,b) = cc * C_orth(a,b)
      dev_d(a,b) = dd * D_orth(a,b)
    enddo ; enddo

    ! Per-corner vertex-based bounds (MLP-u2): corner (a,b) of cell (i,j)
    ! attaches to B-node (I-1+a-1, J-1+b-1) = (i-2+a, j-2+b). Bounds are
    ! relaxed by a small tolerance (Venkatakrishnan 1993 style) so that
    ! machine-precision overshoots at smooth local extrema of the cell-
    ! mean field don't spuriously fire the limiter. The tolerance is
    ! max(BOUND_TOL_ABS, BOUND_TOL_REL * envelope_width) per corner; both
    ! constants are below any physically meaningful overshoot.
    do b = 1, 2 ; do a = 1, 2
      I_node = i - 2 + a ; J_node = j - 2 + b
      if (count_B(I_node, J_node) >= 1.5 .and. &
          Hmax_B(I_node, J_node) > -H_LARGE + 1.0 .and. &
          Hmin_B(I_node, J_node) <  H_LARGE - 1.0) then
        bound_tol = max(BOUND_TOL_ABS, &
                        BOUND_TOL_REL * (Hmax_B(I_node, J_node) - Hmin_B(I_node, J_node)))
        ! Park-Kim MLP-u2 smooth-extrema expansion: relax the envelope by
        ! the theoretical Q1 smooth-peak overshoot magnitude. From the
        ! Taylor expansion of a quadratic field, corner - bound_max at a
        ! smooth peak is ~ |d2|/12. We use (|d2x| + |d2y|)/12 plus a small
        ! fraction of the envelope width as the slack. Cells flagged as
        ! smooth (pk_factor = 1) get the full slack; oscillatory cells
        ! (pk_factor = 0) get strict MLP-u2 (slack = 0).
        pk_slack = pk_factor(i, j) * (PK_SLACK_D2 * (abs(pk_d2x_cell(i, j)) + abs(pk_d2y_cell(i, j))) + &
                                       PK_SLACK_ENV * (Hmax_B(I_node, J_node) - Hmin_B(I_node, J_node)))
        bound_max(a,b) = Hmax_B(I_node, J_node) + bound_tol + pk_slack
        bound_min(a,b) = Hmin_B(I_node, J_node) - bound_tol - pk_slack
      else
        bound_max(a,b) =  H_LARGE
        bound_min(a,b) = -H_LARGE
      endif
    enddo ; enddo

    if (CS%dg_lim_isotropic) then
      ! Isotropic single-phi Zhang-Shu MPP (Zhang & Shu 2010 JCP 229) with
      ! Park-Kim MLP-u2 vertex envelopes. Scales the entire deviation
      ! (h_nodal - Hbar) at each corner by a single per-cell phi, found as
      ! the largest factor in [0,1] that keeps every corner within its
      ! relaxed envelope. Mass conservation is automatic: the deviation
      ! field has zero cell_mean_w-weighted mean by construction.
      phi_iso = 1.0
      do b = 1, 2 ; do a = 1, 2
        iso_dev = CS%h_nodal(i,j,a,b) - Hbar
        if (abs(iso_dev) > TINY_DEV) then
          if (iso_dev > 0.0) then
            ratio = (bound_max(a,b) - Hbar) / iso_dev
          else
            ratio = (bound_min(a,b) - Hbar) / iso_dev
          endif
          phi_iso = min(phi_iso, max(0.0, ratio))
        endif
      enddo ; enddo
      phi_iso = min(phi_iso, 1.0)

      do b = 1, 2 ; do a = 1, 2
        h_new(a,b) = Hbar + phi_iso * (CS%h_nodal(i,j,a,b) - Hbar)
      enddo ; enddo

      CS%h_nodal(i,j,1,1) = h_new(1,1) ; CS%h_nodal(i,j,2,1) = h_new(2,1)
      CS%h_nodal(i,j,1,2) = h_new(1,2) ; CS%h_nodal(i,j,2,2) = h_new(2,2)

      ! Per-mode phi diagnostics all equal in isotropic mode.
      if (associated(CS%dg_lim_phi_xi))    CS%dg_lim_phi_xi(i,j)    = phi_iso
      if (associated(CS%dg_lim_phi_eta))   CS%dg_lim_phi_eta(i,j)   = phi_iso
      if (associated(CS%dg_lim_phi_cross)) CS%dg_lim_phi_cross(i,j) = phi_iso
      if (associated(CS%dg_lim_phi))       CS%dg_lim_phi(i,j)       = phi_iso
      if (associated(CS%dg_lim_mass_drift)) then
        Hbar_new = nodal_cell_mean(h_new, CS%cell_mean_w(i,j,:,:))
        CS%dg_lim_mass_drift(i,j) = Hbar_new - Hbar
      endif
      cycle  ! Skip the anisotropic max-product block below.
    endif

    ! Max-product objective for per-mode scaling. Maximises
    !   phi_b * phi_c * phi_d
    ! (equivalently log(phi_b) + log(phi_c) + log(phi_d)) subject to the 8
    ! corner inequality constraints and the 0 <= phi_i <= 1 box. Unlike a
    ! linear (max-sum) objective, this is symmetric and CONCAVE, so its
    ! optima are graduated (interior of the box where feasible) rather
    ! than vertex-attracted. Each mode is reduced in proportion to its
    ! contribution to bound violation; a mode that contributes nothing
    ! stays at 1, a mode that contributes half of the overshoot gets
    ! half of the reduction.
    !
    ! Solved by iterative single-constraint projection. For one active
    ! upper-bound corner constraint A(:,k) . phi <= B(k), KKT gives:
    !   phi_i = B_eff / (n_pos * A(i,k))  for modes with A(i,k) > 0,
    !   phi_i = 1                          for modes with A(i,k) <= 0,
    ! where B_eff = B(k) - sum of A(i,k) over modes pinned at phi_i = 1.
    ! For multiple binding constraints we iterate: at each step project
    ! onto the most-violated constraint with the KKT formula, clip to
    ! [0, 1], and re-check feasibility. Converges in a few iterations
    ! for typical ice-shelf cells with at most a couple of binding
    ! constraints simultaneously.
    !
    ! Encoding: constraints 1-4 are corner upper bounds, 5-8 are corner
    ! lower bounds (negated). Box constraints are enforced by post-step
    ! clipping rather than as explicit constraints.
    do b = 1, 2 ; do a = 1, 2
      lp_k = (b - 1) * 2 + a
      lp_A(1, lp_k) =  dev_b(a,b)
      lp_A(2, lp_k) =  dev_c(a,b)
      lp_A(3, lp_k) =  dev_d(a,b)
      lp_B(   lp_k) =  bound_max(a,b) - Hbar
      lp_A(1, lp_k+4) = -dev_b(a,b)
      lp_A(2, lp_k+4) = -dev_c(a,b)
      lp_A(3, lp_k+4) = -dev_d(a,b)
      lp_B(   lp_k+4) = -(bound_min(a,b) - Hbar)
    enddo ; enddo

    phi_b = 1.0 ; phi_c = 1.0 ; phi_d = 1.0
    do lp_q = 1, MP_MAX_ITER
      ! Find most-violated constraint at the current phi.
      lp_best_obj = LP_TOL
      lp_i = -1
      do lp_k = 1, 8
        lp_obj = ((lp_A(1,lp_k)*phi_b + lp_A(2,lp_k)*phi_c) + lp_A(3,lp_k)*phi_d) - lp_B(lp_k)
        if (lp_obj > lp_best_obj) then
          lp_best_obj = lp_obj ; lp_i = lp_k
        endif
      enddo
      if (lp_i < 0) exit  ! All constraints satisfied -> feasible.

      ! Single-constraint max-product projection (KKT). Sum positive-coef
      ! contributions and the "ballast" from modes pinned at 1.
      lp_x(1) = phi_b ; lp_x(2) = phi_c ; lp_x(3) = phi_d  ! Save in case of fallback.
      lp_inv_det = lp_B(lp_i)  ! Repurpose lp_inv_det as B_eff working accumulator.
      mp_n_pos = 0
      do lp_k = 1, 3
        if (lp_A(lp_k, lp_i) > 0.0) then
          mp_n_pos = mp_n_pos + 1
        else
          ! Pin this mode at 1; subtract its full contribution from B_eff.
          lp_inv_det = lp_inv_det - lp_A(lp_k, lp_i) * 1.0
          if (lp_k == 1) phi_b = 1.0
          if (lp_k == 2) phi_c = 1.0
          if (lp_k == 3) phi_d = 1.0
        endif
      enddo
      if (mp_n_pos == 0) exit  ! No way to reduce positive-coef modes; nothing helps.

      ! Apply KKT formula to positive-coef modes; clip to [0, 1].
      if (lp_inv_det > 0.0) then
        do lp_k = 1, 3
          if (lp_A(lp_k, lp_i) > 0.0) then
            lp_obj = lp_inv_det / (real(mp_n_pos) * lp_A(lp_k, lp_i))
            lp_obj = max(0.0, min(1.0, lp_obj))
            if (lp_k == 1) phi_b = min(phi_b, lp_obj)
            if (lp_k == 2) phi_c = min(phi_c, lp_obj)
            if (lp_k == 3) phi_d = min(phi_d, lp_obj)
          endif
        enddo
      else
        ! Constraint infeasible at any phi_i >= 0 with positive coefs.
        ! Set positive-coef modes to 0.
        do lp_k = 1, 3
          if (lp_A(lp_k, lp_i) > 0.0) then
            if (lp_k == 1) phi_b = 0.0
            if (lp_k == 2) phi_c = 0.0
            if (lp_k == 3) phi_d = 0.0
          endif
        enddo
      endif
    enddo
    phi_b = max(0.0, min(1.0, phi_b))
    phi_c = max(0.0, min(1.0, phi_c))
    phi_d = max(0.0, min(1.0, phi_d))

    ! Reconstruct corners using orthogonalised modes (preserves Hbar exactly
    ! because each orthogonalised mode has zero cell_mean_w-weighted mean).
    do b = 1, 2 ; do a = 1, 2
      h_new(a,b) = Hbar + phi_b * dev_b(a,b) + phi_c * dev_c(a,b) + phi_d * dev_d(a,b)
    enddo ; enddo

    CS%h_nodal(i,j,1,1) = h_new(1,1) ; CS%h_nodal(i,j,2,1) = h_new(2,1)
    CS%h_nodal(i,j,1,2) = h_new(1,2) ; CS%h_nodal(i,j,2,2) = h_new(2,2)

    ! Diagnostics.
    if (associated(CS%dg_lim_phi_xi))    CS%dg_lim_phi_xi(i,j)    = phi_b
    if (associated(CS%dg_lim_phi_eta))   CS%dg_lim_phi_eta(i,j)   = phi_c
    if (associated(CS%dg_lim_phi_cross)) CS%dg_lim_phi_cross(i,j) = phi_d
    if (associated(CS%dg_lim_phi))       CS%dg_lim_phi(i,j)       = min(phi_b, min(phi_c, phi_d))
    if (associated(CS%dg_lim_mass_drift)) then
      Hbar_new = nodal_cell_mean(h_new, CS%cell_mean_w(i,j,:,:))
      CS%dg_lim_mass_drift(i,j) = Hbar_new - Hbar
    endif
  enddo ; enddo

  call pass_corner_field(CS%h_nodal, G)
end subroutine nodal_hierarchical_limit

!> Subgrid (Phisub) cell-mean of surface elevation s(x,y) = h - bed (grounded) or
!! (1 - rho_i/rho_w)*h (floating), with per-sub-QP flotation determination. Used
!! at grounding-line cells so the cell-mean of s respects partial grounding the
!! same way the strong-form driving stress integration does.
pure real function subgrid_cell_mean_s(Phisub, h_nodal_cell, bed_corners, &
                                       rhoi_rhow, min_h_shelf) result(Sbar)
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub  !< Sub-grid quadrature weights [nondim]
  real, dimension(2,2),         intent(in) :: h_nodal_cell !< Q1 nodal thickness [Z ~> m]
  real, dimension(2,2),         intent(in) :: bed_corners  !< Bed depth at the 4 cell corners
                                                           !! [Z ~> m]
  real, intent(in) :: rhoi_rhow   !< Ice/ocean density ratio [nondim]
  real, intent(in) :: min_h_shelf !< Positivity floor for thickness [Z ~> m]

  real :: h_gp        ! Bilinear-interpolated thickness at a sub-QP [Z ~> m]
  real :: bed_gp      ! Bilinear-interpolated bed depth at a sub-QP [Z ~> m]
  real :: s_gp        ! Surface elevation at a sub-QP [Z ~> m]
  real :: subarea     ! Reference-cell area of one sub-cell [nondim]
  real :: accum       ! Running area-weighted sum of s [Z ~> m]
  integer :: nsub, ii, jj, qx, qy
  logical :: is_grounded

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  ! sum over all (ii,jj,qx,qy) of (0.25 * subarea * s_gp) gives the area-weighted
  ! cell-mean of s on the reference cell [0,1]^2 (sum of weights = 1).
  accum = 0.0
  do jj = 1, nsub ; do ii = 1, nsub ; do qy = 1, 2 ; do qx = 1, 2
    h_gp = ((Phisub(qx,qy,ii,jj,1,1) * h_nodal_cell(1,1)) + &
            (Phisub(qx,qy,ii,jj,2,2) * h_nodal_cell(2,2))) + &
           ((Phisub(qx,qy,ii,jj,1,2) * h_nodal_cell(1,2)) + &
            (Phisub(qx,qy,ii,jj,2,1) * h_nodal_cell(2,1)))
    h_gp = max(h_gp, min_h_shelf)
    bed_gp = ((Phisub(qx,qy,ii,jj,1,1) * bed_corners(1,1)) + &
              (Phisub(qx,qy,ii,jj,2,2) * bed_corners(2,2))) + &
             ((Phisub(qx,qy,ii,jj,1,2) * bed_corners(1,2)) + &
              (Phisub(qx,qy,ii,jj,2,1) * bed_corners(2,1)))
    is_grounded = (rhoi_rhow * h_gp - bed_gp > 0.0)
    if (is_grounded) then
      s_gp = h_gp - bed_gp
    else
      s_gp = (1.0 - rhoi_rhow) * h_gp
    endif
    accum = accum + (0.25 * subarea) * s_gp
  enddo ; enddo ; enddo ; enddo
  Sbar = accum
end function subgrid_cell_mean_s

!> Well-balanced surface-elevation limiter for DG(1) ice-shelf thickness.
!! Limits the surface elevation s instead of thickness h, matching the
!! well-balancing principle from shallow-water DG (Audusse 2004 hydrostatic
!! reconstruction, Xing-Shu 2005 well-balanced DG): on rough bedrock the
!! equilibrium state has dh/dx ~ dbed/dx, so limiting h destroys the
!! compensation and produces spurious surface slopes equal to the bed slope.
!! Limiting s preserves the dynamically smooth variable directly.
!!
!! Per-corner surface elevation:
!!   grounded corner (rhoi_rhow*h > bed):  s = h - bed
!!   floating corner:                       s = (1 - rhoi_rhow)*h
!! Flotation guarantees the formulas agree at the GL so s is continuous within
!! a cell with mixed corners.
!!
!! Cell-mean of s for the envelope:
!!   fully grounded or fully floating cell: closed-form linear in h.
!!   GL cell (0 < ground_frac < 1): Phisub subgrid quadrature (the same Phisub
!!   used by the strong-form driving stress integration).
!!
!! Algorithm per cell:
!!   1) Build s_nodal at each corner using pre-limit flotation state.
!!   2) Isotropic Zhang-Shu single-phi MPP with Park-Kim MLP-u2 slack on s.
!!   3) Reconstruct h corner-by-corner from limited s using the *frozen*
!!      pre-limit flotation state (avoids mid-step flotation flips).
!!   4) Mass-fix uniform shift: h_final = h_recon + (Hbar_old - Hbar_new) so
!!      cell_mean(h) is restored exactly to its pre-limit value at every
!!      limited cell.
!!
!! References: Audusse et al. 2004 SIAM JSC 25, Xing & Shu 2005 JCP 208,
!! Zhang & Shu 2010 JCP 229, Park & Kim 2014 JCP 274.
subroutine nodal_surface_slope_limit(CS, G, ISS)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ocean_grid_type),  intent(in)    :: G
  type(ice_shelf_state),  intent(in)    :: ISS

  real, dimension(SZDIB_(G),SZDJB_(G)) :: Smax_B, Smin_B   ! Per-B-node max/min of cell-mean
                                                           ! s over touching cells [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: count_B          ! Number of touching valid cells
                                                           ! at each B-node [nondim]
  real, dimension(SZDI_(G),SZDJ_(G))   :: S_cell           ! Cell-mean of s [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G))   :: pk_factor        ! Park-Kim smooth-extrema flag,
                                                           ! 0 or 1 [nondim]
  real, dimension(SZDI_(G),SZDJ_(G))   :: pk_d2x_cell      ! Second x-difference of S_cell
                                                           ! [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G))   :: pk_d2y_cell      ! Second y-difference of S_cell
                                                           ! [Z ~> m]
  real :: pk_slack    ! Park-Kim envelope relaxation [Z ~> m]
  real :: bound_tol   ! Per-corner Venkatakrishnan-style envelope tolerance [Z ~> m]
  real :: env_width   ! Smax_B - Smin_B at a B-node [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G))   :: venkat_slack    ! Per-cell gradient-proportional
                                                          ! slack from centred differences
                                                          ! of S_cell [Z ~> m]
  real :: grad_S_x    ! Centred-difference estimate of cell-mean dS/dx [Z ~> m]
  real :: grad_S_y    ! Centred-difference estimate of cell-mean dS/dy [Z ~> m]
  real :: phi_iso     ! Isotropic per-cell MPP scaling factor [nondim]
  real :: s_dev_ab    ! Per-corner deviation s_nodal - Sbar [Z ~> m]
  real :: ratio       ! Per-corner s-bound / s-deviation ratio [nondim]
  real, parameter :: H_LARGE = 1.0e30      ! Sentinel for "no envelope" [Z ~> m]
  real, parameter :: TINY_DEV = 1.0e-30    ! Skip threshold for near-zero deviation [Z ~> m]
  real, parameter :: PK_SLACK_D2  = 1.0/12.0  ! Q1 smooth-peak overshoot coefficient [nondim]
  real, parameter :: PK_SLACK_ENV = 0.05      ! Envelope-width slack fraction [nondim]
  real, parameter :: BOUND_TOL_ABS = 1.0e-4   ! Absolute envelope floor [Z ~> m]
  real, parameter :: BOUND_TOL_REL = 1.0e-6   ! Relative envelope floor [nondim]
  real :: rhoi_rhow         ! Ice/ocean density ratio rho_i/rho_w [nondim]
  real :: one_minus_r       ! 1 - rhoi_rhow [nondim]
  real :: inv_one_minus_r   ! 1 / (1 - rhoi_rhow) [nondim]
  real :: Sbar              ! Cell-mean surface elevation [Z ~> m]
  real :: Hbar_old          ! Pre-limit cell-mean thickness [Z ~> m]
  real :: Hbar_new          ! Reconstructed cell-mean thickness before mass-fix [Z ~> m]
  real :: shift             ! Mass-fix uniform shift Hbar_old - Hbar_new [Z ~> m]
  real :: bed_c(2,2)        ! Bed depth at the 4 cell corners [Z ~> m]
  logical :: cg(2,2)        ! Per-corner pre-limit grounded flag
  logical :: cell_is_GL     ! True if the cell has mixed grounded/floating corners
  real :: s_nodal(2,2)      ! Per-corner pre-limit surface elevation [Z ~> m]
  real :: s_limited(2,2)    ! Per-corner post-limit surface elevation [Z ~> m]
  real :: s_bound_max(2,2)  ! Per-corner relaxed upper s-bound [Z ~> m]
  real :: s_bound_min(2,2)  ! Per-corner relaxed lower s-bound [Z ~> m]
  real :: h_recon(2,2)      ! Per-corner reconstructed thickness before mass-fix [Z ~> m]
  real :: hbdry_val         ! Dirichlet h_bdry_val floored by min_h_shelf [Z ~> m]
  real :: bed_avg_bdry      ! Average bed depth over 4 corners (Dirichlet cell) [Z ~> m]
  integer :: i, j, a, b, I_node, J_node

  if (.not. associated(CS%h_nodal)) return
  if (.not. associated(CS%Phisub)) return  ! Subgrid weights required for GL cells

  rhoi_rhow       = CS%density_ice / CS%density_ocean_avg
  one_minus_r     = 1.0 - rhoi_rhow
  inv_one_minus_r = 1.0 / one_minus_r

  call pass_corner_field(CS%h_nodal, G)

  ! 1) Build S_cell = cell-mean of surface elevation s over the data domain.
  !    For each cell, classify the cell as fully grounded / fully floating / GL
  !    by the per-corner flotation test on (h_nodal, bed_node). GL cells use
  !    Phisub subgrid quadrature to integrate s respecting subgrid GL position.
  S_cell(:,:) = H_LARGE
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    if (ISS%hmask(i,j) == 1.0) then
      ! Gather bed at corners and per-corner pre-limit flotation state.
      bed_c(1,1) = CS%bed_node(i-1,j-1) ; bed_c(2,1) = CS%bed_node(i,j-1)
      bed_c(1,2) = CS%bed_node(i-1,j  ) ; bed_c(2,2) = CS%bed_node(i,j  )
      cg(1,1) = (rhoi_rhow * CS%h_nodal(i,j,1,1) > bed_c(1,1))
      cg(2,1) = (rhoi_rhow * CS%h_nodal(i,j,2,1) > bed_c(2,1))
      cg(1,2) = (rhoi_rhow * CS%h_nodal(i,j,1,2) > bed_c(1,2))
      cg(2,2) = (rhoi_rhow * CS%h_nodal(i,j,2,2) > bed_c(2,2))
      cell_is_GL = .not. ((cg(1,1) .and. cg(2,1) .and. cg(1,2) .and. cg(2,2)) .or. &
                          (.not.(cg(1,1) .or. cg(2,1) .or. cg(1,2) .or. cg(2,2))))
      if (cell_is_GL) then
        S_cell(i,j) = subgrid_cell_mean_s(CS%Phisub, CS%h_nodal(i,j,:,:), bed_c, &
                                          rhoi_rhow, CS%min_h_shelf)
      else
        do b = 1, 2 ; do a = 1, 2
          if (cg(a,b)) then
            s_nodal(a,b) = CS%h_nodal(i,j,a,b) - bed_c(a,b)
          else
            s_nodal(a,b) = one_minus_r * CS%h_nodal(i,j,a,b)
          endif
        enddo ; enddo
        S_cell(i,j) = nodal_cell_mean(s_nodal, CS%cell_mean_w(i,j,:,:))
      endif
    elseif (ISS%hmask(i,j) == 3.0) then
      ! Dirichlet thickness BC: derive s from h_bdry_val and the average bed at
      ! this cell's 4 corners. Use the same flotation test on the cell-average
      ! bed to choose the formula.
      hbdry_val = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
      bed_avg_bdry = 0.25 * ((CS%bed_node(i-1,j-1) + CS%bed_node(i,j  )) + &
                             (CS%bed_node(i,  j-1) + CS%bed_node(i-1,j)))
      if (rhoi_rhow * hbdry_val > bed_avg_bdry) then
        S_cell(i,j) = hbdry_val - bed_avg_bdry
      else
        S_cell(i,j) = one_minus_r * hbdry_val
      endif
    else
      cycle
    endif
  enddo ; enddo

  ! 2) Scatter S_cell to B-nodes -> Smax_B, Smin_B (MLP-u2 vertex envelope on s).
  Smax_B(:,:) = -H_LARGE
  Smin_B(:,:) =  H_LARGE
  count_B(:,:) = 0.0
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    if (S_cell(i,j) >= H_LARGE - 1.0) cycle
    if (i-1 >= G%IsdB .and. j-1 >= G%JsdB) then
      Smax_B(i-1, j-1) = max(Smax_B(i-1, j-1), S_cell(i,j))
      Smin_B(i-1, j-1) = min(Smin_B(i-1, j-1), S_cell(i,j))
      count_B(i-1, j-1) = count_B(i-1, j-1) + 1.0
    endif
    if (i   <= G%IedB .and. j-1 >= G%JsdB) then
      Smax_B(i,   j-1) = max(Smax_B(i,   j-1), S_cell(i,j))
      Smin_B(i,   j-1) = min(Smin_B(i,   j-1), S_cell(i,j))
      count_B(i,   j-1) = count_B(i,   j-1) + 1.0
    endif
    if (i-1 >= G%IsdB .and. j   <= G%JedB) then
      Smax_B(i-1, j  ) = max(Smax_B(i-1, j  ), S_cell(i,j))
      Smin_B(i-1, j  ) = min(Smin_B(i-1, j  ), S_cell(i,j))
      count_B(i-1, j  ) = count_B(i-1, j  ) + 1.0
    endif
    if (i   <= G%IedB .and. j   <= G%JedB) then
      Smax_B(i,   j  ) = max(Smax_B(i,   j  ), S_cell(i,j))
      Smin_B(i,   j  ) = min(Smin_B(i,   j  ), S_cell(i,j))
      count_B(i,   j  ) = count_B(i,   j  ) + 1.0
    endif
  enddo ; enddo
  call pass_var(Smax_B,  G%domain, position=CORNER)
  call pass_var(Smin_B,  G%domain, position=CORNER)
  call pass_var(count_B, G%domain, position=CORNER)

  ! 3) Park-Kim MLP-u2 smooth-extrema indicator on S_cell. Same 3-cell stencil
  !    sign-consistency check (Park & Kim 2014 JCP 274) but applied to the
  !    surface-elevation field, which is the dynamically smooth variable.
  call pass_var(S_cell, G%domain)
  pk_d2x_cell(:,:) = 0.0
  pk_d2y_cell(:,:) = 0.0
  do j = G%jsd+1, G%jed-1 ; do i = G%isd+1, G%ied-1
    if (S_cell(i, j) >= H_LARGE - 1.0) cycle
    if (S_cell(i-1, j) < H_LARGE - 1.0 .and. S_cell(i+1, j) < H_LARGE - 1.0) then
      pk_d2x_cell(i, j) = (S_cell(i-1, j) - 2.0 * S_cell(i, j)) + S_cell(i+1, j)
    endif
    if (S_cell(i, j-1) < H_LARGE - 1.0 .and. S_cell(i, j+1) < H_LARGE - 1.0) then
      pk_d2y_cell(i, j) = (S_cell(i, j-1) - 2.0 * S_cell(i, j)) + S_cell(i, j+1)
    endif
  enddo ; enddo
  ! Venkatakrishnan-style per-cell gradient-proportional slack from centred
  ! differences of S_cell. Corner-vs-Sbar offset of a linear Q1 field with
  ! cell-mean gradients (grad_S_x, grad_S_y) is (|grad_S_x|+|grad_S_y|)/2
  ! when grad_S is expressed per cell width. We estimate grad per cell by the
  ! centred difference of neighbour cell means (= 2 cells apart), giving
  ! grad_S_x ~ (S(i+1)-S(i-1))/2 in per-cell units; halved again for the
  ! half-cell-width corner offset gives /4. K = DG1_VENKAT_K scales the term;
  ! K = 1 grants exactly the linear-field corner offset.
  venkat_slack(:,:) = 0.0
  do j = G%jsd+1, G%jed-1 ; do i = G%isd+1, G%ied-1
    if (S_cell(i, j) >= H_LARGE - 1.0) cycle
    grad_S_x = 0.0
    grad_S_y = 0.0
    if (S_cell(i-1, j) < H_LARGE - 1.0 .and. S_cell(i+1, j) < H_LARGE - 1.0) &
      grad_S_x = S_cell(i+1, j) - S_cell(i-1, j)
    if (S_cell(i, j-1) < H_LARGE - 1.0 .and. S_cell(i, j+1) < H_LARGE - 1.0) &
      grad_S_y = S_cell(i, j+1) - S_cell(i, j-1)
    venkat_slack(i, j) = CS%dg_venkat_K * 0.25 * (abs(grad_S_x) + abs(grad_S_y))
  enddo ; enddo

  pk_factor(:,:) = 0.0
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle
    if (i-1 < G%isd+1 .or. i+1 > G%ied-1) cycle
    if (j-1 < G%jsd+1 .or. j+1 > G%jed-1) cycle
    if (((pk_d2x_cell(i-1, j) >= 0.0) .and. (pk_d2x_cell(i, j) >= 0.0) .and. &
         (pk_d2x_cell(i+1, j) >= 0.0)) .or. &
        ((pk_d2x_cell(i-1, j) <= 0.0) .and. (pk_d2x_cell(i, j) <= 0.0) .and. &
         (pk_d2x_cell(i+1, j) <= 0.0))) then
      if (((pk_d2y_cell(i, j-1) >= 0.0) .and. (pk_d2y_cell(i, j) >= 0.0) .and. &
           (pk_d2y_cell(i, j+1) >= 0.0)) .or. &
          ((pk_d2y_cell(i, j-1) <= 0.0) .and. (pk_d2y_cell(i, j) <= 0.0) .and. &
           (pk_d2y_cell(i, j+1) <= 0.0))) then
        pk_factor(i, j) = 1.0
      endif
    endif
  enddo ; enddo

  ! Reset diagnostic buffers (limiter inactive at non-hmask=1 cells -> phi=1, drift=0).
  if (associated(CS%dg_lim_phi_xi))     CS%dg_lim_phi_xi(:,:)     = 1.0
  if (associated(CS%dg_lim_phi_eta))    CS%dg_lim_phi_eta(:,:)    = 1.0
  if (associated(CS%dg_lim_phi_cross))  CS%dg_lim_phi_cross(:,:)  = 1.0
  if (associated(CS%dg_lim_mass_drift)) CS%dg_lim_mass_drift(:,:) = 0.0
  if (associated(CS%dg_lim_phi))        CS%dg_lim_phi(:,:)        = 1.0
  if (associated(CS%dg_lim_pk_factor))  CS%dg_lim_pk_factor(:,:)  = pk_factor(:,:)

  ! 4) Per-cell isotropic limit on s, reconstruct h, mass-fix uniform shift.
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle

    ! Snapshot pre-limit h-cell-mean for the mass-fix shift below.
    Hbar_old = nodal_cell_mean(CS%h_nodal(i,j,:,:), CS%cell_mean_w(i,j,:,:))

    ! Bed at corners and frozen per-corner flotation state.
    bed_c(1,1) = CS%bed_node(i-1,j-1) ; bed_c(2,1) = CS%bed_node(i,j-1)
    bed_c(1,2) = CS%bed_node(i-1,j  ) ; bed_c(2,2) = CS%bed_node(i,j  )
    cg(1,1) = (rhoi_rhow * CS%h_nodal(i,j,1,1) > bed_c(1,1))
    cg(2,1) = (rhoi_rhow * CS%h_nodal(i,j,2,1) > bed_c(2,1))
    cg(1,2) = (rhoi_rhow * CS%h_nodal(i,j,1,2) > bed_c(1,2))
    cg(2,2) = (rhoi_rhow * CS%h_nodal(i,j,2,2) > bed_c(2,2))
    cell_is_GL = .not. ((cg(1,1) .and. cg(2,1) .and. cg(1,2) .and. cg(2,2)) .or. &
                        (.not.(cg(1,1) .or. cg(2,1) .or. cg(1,2) .or. cg(2,2))))

    ! s_nodal at each corner.
    do b = 1, 2 ; do a = 1, 2
      if (cg(a,b)) then
        s_nodal(a,b) = CS%h_nodal(i,j,a,b) - bed_c(a,b)
      else
        s_nodal(a,b) = one_minus_r * CS%h_nodal(i,j,a,b)
      endif
    enddo ; enddo

    Sbar = S_cell(i,j)

    ! Per-corner relaxed s-envelope (MLP-u2 + Park-Kim smooth-extrema slack).
    do b = 1, 2 ; do a = 1, 2
      I_node = i - 2 + a ; J_node = j - 2 + b
      if (count_B(I_node, J_node) >= 1.5 .and. &
          Smax_B(I_node, J_node) > -H_LARGE + 1.0 .and. &
          Smin_B(I_node, J_node) <  H_LARGE - 1.0) then
        env_width = Smax_B(I_node, J_node) - Smin_B(I_node, J_node)
        bound_tol = max(BOUND_TOL_ABS, BOUND_TOL_REL * env_width)
        pk_slack = pk_factor(i, j) * &
                   (PK_SLACK_D2 * (abs(pk_d2x_cell(i, j)) + abs(pk_d2y_cell(i, j))) + &
                    PK_SLACK_ENV * env_width)
        s_bound_max(a,b) = Smax_B(I_node, J_node) + bound_tol + pk_slack + venkat_slack(i, j)
        s_bound_min(a,b) = Smin_B(I_node, J_node) - bound_tol - pk_slack - venkat_slack(i, j)
      else
        s_bound_max(a,b) =  H_LARGE
        s_bound_min(a,b) = -H_LARGE
      endif
    enddo ; enddo

    ! Isotropic single-phi Zhang-Shu MPP on s.
    phi_iso = 1.0
    do b = 1, 2 ; do a = 1, 2
      s_dev_ab = s_nodal(a,b) - Sbar
      if (abs(s_dev_ab) > TINY_DEV) then
        if (s_dev_ab > 0.0) then
          ratio = (s_bound_max(a,b) - Sbar) / s_dev_ab
        else
          ratio = (s_bound_min(a,b) - Sbar) / s_dev_ab
        endif
        phi_iso = min(phi_iso, max(0.0, ratio))
      endif
    enddo ; enddo
    phi_iso = min(phi_iso, 1.0)

    if (cell_is_GL) then
      ! GL cells: apply phi (computed from s-envelope) directly to h to avoid the
      ! s<->h Jacobian discontinuity at the in-cell GL. Mass is trivially
      ! conserved (linear scaling of deviation from Hbar). Well-balanced
      ! reconstruction is given up inside GL cells only; the s-based firing
      ! decision still benefits from the smooth-surface criterion.
      do b = 1, 2 ; do a = 1, 2
        h_recon(a,b) = Hbar_old + phi_iso * (CS%h_nodal(i,j,a,b) - Hbar_old)
      enddo ; enddo
      if (associated(CS%dg_lim_mass_drift)) CS%dg_lim_mass_drift(i,j) = 0.0
      CS%h_nodal(i,j,1,1) = h_recon(1,1) ; CS%h_nodal(i,j,2,1) = h_recon(2,1)
      CS%h_nodal(i,j,1,2) = h_recon(1,2) ; CS%h_nodal(i,j,2,2) = h_recon(2,2)
    else
      ! Non-GL cells: full well-balanced reconstruction. Limit s, recover h
      ! corner-by-corner from the frozen flotation state, then mass-fix uniform
      ! shift to restore cell_mean(h) exactly (handles FP roundoff in
      ! nodal_cell_mean; the algorithmic drift is zero by linearity on
      ! uniform-formula cells).
      do b = 1, 2 ; do a = 1, 2
        s_limited(a,b) = Sbar + phi_iso * (s_nodal(a,b) - Sbar)
        if (cg(a,b)) then
          h_recon(a,b) = s_limited(a,b) + bed_c(a,b)
        else
          h_recon(a,b) = s_limited(a,b) * inv_one_minus_r
        endif
      enddo ; enddo
      Hbar_new = nodal_cell_mean(h_recon, CS%cell_mean_w(i,j,:,:))
      if (associated(CS%dg_lim_mass_drift)) CS%dg_lim_mass_drift(i,j) = Hbar_new - Hbar_old
      shift = Hbar_old - Hbar_new
      CS%h_nodal(i,j,1,1) = h_recon(1,1) + shift
      CS%h_nodal(i,j,2,1) = h_recon(2,1) + shift
      CS%h_nodal(i,j,1,2) = h_recon(1,2) + shift
      CS%h_nodal(i,j,2,2) = h_recon(2,2) + shift
    endif

    if (associated(CS%dg_lim_phi_xi))    CS%dg_lim_phi_xi(i,j)    = phi_iso
    if (associated(CS%dg_lim_phi_eta))   CS%dg_lim_phi_eta(i,j)   = phi_iso
    if (associated(CS%dg_lim_phi_cross)) CS%dg_lim_phi_cross(i,j) = phi_iso
    if (associated(CS%dg_lim_phi))       CS%dg_lim_phi(i,j)       = phi_iso
  enddo ; enddo

  call pass_corner_field(CS%h_nodal, G)
end subroutine nodal_surface_slope_limit

!> Apply the per-cell tensor-product Q1 mass-matrix inverse:
!! out(a,b) = sum_{a',b'} Minv_xi(a,a') * Minv_eta(b,b') * rhs(a',b').
pure subroutine apply_nodal_DG_mass_inverse(Minv_xi_cell, Minv_eta_cell, rhs, out)
  real, dimension(2,2), intent(in)  :: Minv_xi_cell  !< 2x2 inverse of xi mass-matrix factor [L-1]
  real, dimension(2,2), intent(in)  :: Minv_eta_cell !< 2x2 inverse of eta mass-matrix factor [L-1]
  real, dimension(2,2), intent(in)  :: rhs            !< Per-cell RHS at the 4 corners [Z L2 T-1]
  real, dimension(2,2), intent(out) :: out            !< M^-1 * rhs [Z T-1]
  integer :: a, b, ap, bp
  real :: s
  do b = 1, 2 ; do a = 1, 2
    s = 0.0
    do bp = 1, 2 ; do ap = 1, 2
      s = s + Minv_xi_cell(a,ap) * Minv_eta_cell(b,bp) * rhs(ap,bp)
    enddo ; enddo
    out(a,b) = s
  enddo ; enddo
end subroutine apply_nodal_DG_mass_inverse

!> Per-side-flotation surface-elevation jump s_B - s_A at a face quadrature point.
!! Each side evaluates its own flotation branch on its own thickness and the shared
!! single-valued face bed. A hydrostatically-continuous grounding line (s_B = s_A)
!! returns zero even when h_B /= h_A.
pure function dg1_wb_surface_jump(h_A, h_B, bed_qp, rhoi_rhow) result(ds)
  real, intent(in) :: h_A       !< Side-A face-QP thickness [Z ~> m]
  real, intent(in) :: h_B       !< Side-B face-QP thickness [Z ~> m]
  real, intent(in) :: bed_qp    !< Bed elevation at the face QP, single-valued [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio [nondim]
  real :: ds                    !< Surface-elevation jump s_B - s_A [Z ~> m]
  real :: s_A, s_B              ! Per-side surface elevation [Z ~> m]
  real :: one_m_r               ! 1 - rhoi_rhow [nondim]
  one_m_r = 1.0 - rhoi_rhow
  if (rhoi_rhow*h_A - bed_qp > 0.0) then ; s_A = h_A - bed_qp
  else ; s_A = one_m_r*h_A ; endif
  if (rhoi_rhow*h_B - bed_qp > 0.0) then ; s_B = h_B - bed_qp
  else ; s_B = one_m_r*h_B ; endif
  ds = s_B - s_A
end function dg1_wb_surface_jump

!> Mean-supported surface-jump allowance for the DG(1) excess-jump artificial
!! viscosity: the largest |s_B - s_A| the two cell means can support at the face-QP
!! bed over the admissible flotation-branch assignments. A cell mean's branch test
!! at the face bed can misclassify when the face bed is a local extremum
!! unrepresentative of the cell interior (e.g. a grounding line on a sill, where a
!! floating cell's mean reads grounded at the shallow face bed and the allowance
!! collapses to a plain thickness contrast). The admissible branch set per side is
!! {branch of the mean, branch of the face trace}, both tested at the face bed, so
!! the allowance is unchanged wherever mean and trace agree on the branch and only
!! widens (conservatively, toward under-damping) at flotation-ambiguous faces.
pure function dg1_wb_mean_allowance(Hbar_A, Hbar_B, h_A, h_B, bed_qp, rhoi_rhow, &
                                    branch_max) result(ds_allow)
  real, intent(in) :: Hbar_A    !< Side-A cell-mean thickness [Z ~> m]
  real, intent(in) :: Hbar_B    !< Side-B cell-mean thickness [Z ~> m]
  real, intent(in) :: h_A       !< Side-A face-QP trace thickness [Z ~> m]
  real, intent(in) :: h_B       !< Side-B face-QP trace thickness [Z ~> m]
  real, intent(in) :: bed_qp    !< Bed elevation at the face QP, single-valued [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio [nondim]
  logical, intent(in) :: branch_max !< If true take the max over admissible branch
                                !! assignments; if false use the mean's own branch
                                !! only (legacy)
  real :: ds_allow              !< Nonnegative mean-supported allowance on |[s]| [Z ~> m]
  real :: s_A(2), s_B(2)        ! Per-side candidate mean surfaces [Z ~> m]
  real :: one_m_r               ! 1 - rhoi_rhow [nondim]
  integer :: nA, nB, ia, ib
  logical :: gA_mean, gB_mean, gA_trace, gB_trace

  if (.not. branch_max) then
    ds_allow = abs(dg1_wb_surface_jump(Hbar_A, Hbar_B, bed_qp, rhoi_rhow))
    return
  endif

  one_m_r = 1.0 - rhoi_rhow
  gA_mean  = (rhoi_rhow*Hbar_A - bed_qp > 0.0)
  gB_mean  = (rhoi_rhow*Hbar_B - bed_qp > 0.0)
  gA_trace = (rhoi_rhow*h_A - bed_qp > 0.0)
  gB_trace = (rhoi_rhow*h_B - bed_qp > 0.0)

  if (gA_mean) then ; s_A(1) = Hbar_A - bed_qp ; else ; s_A(1) = one_m_r*Hbar_A ; endif
  nA = 1
  if (gA_trace .neqv. gA_mean) then
    nA = 2
    if (gA_trace) then ; s_A(2) = Hbar_A - bed_qp ; else ; s_A(2) = one_m_r*Hbar_A ; endif
  endif
  if (gB_mean) then ; s_B(1) = Hbar_B - bed_qp ; else ; s_B(1) = one_m_r*Hbar_B ; endif
  nB = 1
  if (gB_trace .neqv. gB_mean) then
    nB = 2
    if (gB_trace) then ; s_B(2) = Hbar_B - bed_qp ; else ; s_B(2) = one_m_r*Hbar_B ; endif
  endif

  ds_allow = 0.0
  do ib = 1, nB ; do ia = 1, nA
    ds_allow = max(ds_allow, abs(s_B(ib) - s_A(ia)))
  enddo ; enddo
end function dg1_wb_mean_allowance

!> Mean inverse flotation slope dh/ds across a face (1 grounded, 1/(1-r) floating per
!! side), harmonic or arithmetic. Multiplying a surface jump by this mean gives the
!! equivalent thickness jump. The same-branch shortcut keeps uniform-flotation faces
!! bitwise independent of the mean choice (both means coincide there). The harmonic
!! mean makes the jump-mode stability rate amplification exactly 2 at every face; the
!! arithmetic mean (legacy) amplifies mixed grounded/floating faces by up to ~2.8x
!! relative to that (rate factor 0.5*(1+1/(1-r))*(1+(1-r)) ~ 5.7 vs 2, for r ~= 0.89).
pure function dg1_wb_slope_mean(h_A, h_B, bed_qp, rhoi_rhow, harmonic) result(m)
  real, intent(in) :: h_A       !< Side-A face-QP thickness [Z ~> m]
  real, intent(in) :: h_B       !< Side-B face-QP thickness [Z ~> m]
  real, intent(in) :: bed_qp    !< Bed elevation at the face QP, single-valued [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio [nondim]
  logical, intent(in) :: harmonic !< If true use the harmonic mean of per-side dh/ds,
                                !! else the arithmetic mean
  real :: m                     !< Mean inverse flotation slope dh/ds [nondim]
  real :: g_A, g_B              ! Per-side surface slope ds/dh [nondim]
  real :: one_m_r               ! 1 - rhoi_rhow [nondim]
  one_m_r = 1.0 - rhoi_rhow
  if (rhoi_rhow*h_A - bed_qp > 0.0) then ; g_A = 1.0 ; else ; g_A = one_m_r ; endif
  if (rhoi_rhow*h_B - bed_qp > 0.0) then ; g_B = 1.0 ; else ; g_B = one_m_r ; endif
  if (g_A == g_B) then
    m = 1.0/g_A   ! Same flotation branch: both means coincide; this form keeps
                  ! uniform-flotation faces bitwise identical under either choice.
  elseif (harmonic) then
    m = 2.0/(g_A + g_B)
  else
    m = 0.5*((1.0/g_A) + (1.0/g_B))
  endif
end function dg1_wb_slope_mean

!> Well-balanced equivalent thickness jump for the DG(1) artificial viscosity: the
!! per-side surface-elevation jump [s] mapped back to a thickness flux by the mean
!! inverse flotation slope. On a uniform-flotation face the single-valued bed makes
!! [s] proportional to [h], so the result equals (h_B - h_A) exactly; the harmonic
!! and arithmetic means diverge only at mixed-flotation (grounding-line) faces, where
!! a hydrostatically-continuous surface (s_B = s_A) returns zero even when h_B /= h_A.
pure function dg1_wb_equiv_jump(h_A, h_B, bed_qp, rhoi_rhow, harmonic) result(dh_eq)
  real, intent(in) :: h_A       !< Side-A face-QP thickness [Z ~> m]
  real, intent(in) :: h_B       !< Side-B face-QP thickness [Z ~> m]
  real, intent(in) :: bed_qp    !< Bed elevation at the face QP, single-valued [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio [nondim]
  logical, intent(in) :: harmonic !< If true use the harmonic mean of per-side dh/ds,
                                !! else the arithmetic mean
  real :: dh_eq                 !< Well-balanced equivalent thickness jump [Z ~> m]
  dh_eq = dg1_wb_surface_jump(h_A, h_B, bed_qp, rhoi_rhow) * &
          dg1_wb_slope_mean(h_A, h_B, bed_qp, rhoi_rhow, harmonic)
end function dg1_wb_equiv_jump

!> Jump-mode rate amplification factor for the per-face stability budget of the
!! DG(1) artificial viscosity. The semi-discrete decay rate of the face-node jump
!! mode is lambda = 4 * amp * c * u_eff / dx_perp, where the factor 4 collects the
!! node-localization (x2) and consistent-mass (x2) amplifications relative to the
!! cell-mean rate, and amp = 0.5*(dh/ds|_A + dh/ds|_B) * (ds/dh|_A + ds/dh|_B)
!! collects the flotation-branch flux Jacobian. With the harmonic slope mean amp = 2
!! exactly at every face by construction. With the arithmetic mean amp = 2 on
!! uniform-flotation faces and rises to ~5.7 (for r ~= 0.89) at mixed
!! grounded/floating faces, where the well-balanced flux responds more strongly to a
!! thickness perturbation than the raw jump [h] suggests.
pure function dg1_wb_jump_rate_amp(h_A, h_B, bed_qp, rhoi_rhow, harmonic) result(amp)
  real, intent(in) :: h_A       !< Side-A face-QP thickness [Z ~> m]
  real, intent(in) :: h_B       !< Side-B face-QP thickness [Z ~> m]
  real, intent(in) :: bed_qp    !< Bed elevation at the face QP, single-valued [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio [nondim]
  logical, intent(in) :: harmonic !< If true the harmonic slope mean is in use
  real :: amp                   !< Rate amplification factor, >= 2 [nondim]
  real :: g_A, g_B              ! Per-side surface slope ds/dh [nondim]
  real :: one_m_r               ! 1 - rhoi_rhow [nondim]
  if (harmonic) then
    amp = 2.0   ! (2/(g_A+g_B)) * (g_A+g_B): exact at every face by construction.
  else
    one_m_r = 1.0 - rhoi_rhow
    if (rhoi_rhow*h_A - bed_qp > 0.0) then ; g_A = 1.0 ; else ; g_A = one_m_r ; endif
    if (rhoi_rhow*h_B - bed_qp > 0.0) then ; g_B = 1.0 ; else ; g_B = one_m_r ; endif
    amp = 0.5 * ((1.0/g_A) + (1.0/g_B)) * (g_A + g_B)
  endif
end function dg1_wb_jump_rate_amp

!> 2D effective strain rate (second invariant) used by the DG(1) artificial-viscosity
!! strain-scaled floor. eps_e = sqrt(eps_xx^2 + eps_yy^2 + eps_xx*eps_yy + eps_xy^2) with
!! eps_xy = 0.5*(du/dy + dv/dx). max(0,.) guards FP roundoff in the radicand.
pure function dg1_face_eps_eff(dudx, dudy, dvdx, dvdy) result(eps_e)
  real, intent(in) :: dudx     !< du/dx at face midpoint [T-1]
  real, intent(in) :: dudy     !< du/dy at face midpoint [T-1]
  real, intent(in) :: dvdx     !< dv/dx at face midpoint [T-1]
  real, intent(in) :: dvdy     !< dv/dy at face midpoint [T-1]
  real :: eps_e                !< Effective strain rate [T-1]
  real :: eps_xy               ! Off-diagonal symmetric strain rate [T-1]
  eps_xy = 0.5*(dudy + dvdx)
  eps_e = sqrt(max(0.0, dudx*dudx + dvdy*dvdy + dudx*dvdy + eps_xy*eps_xy))
end function dg1_face_eps_eff

!> Compute the nodal Q1 DG(1) spatial operator (RHS of the per-cell mass-matrix
!! system for d h_nodal / dt). Implements the IBP weak form
!!   M dh/dt = + int grad(N) . (u h) dV  -  contour N (u.n) h_upwind ds
!! at 2x2 Gauss-Legendre QPs in the volume and 2-point Gauss on each face.
!! Specified-flux faces (u/v_face_mask == 4) distribute the prescribed face
!! flux to the two on-face corners with equal weight (plan R24).
subroutine DG1_nodal_spatial_operator(CS, G, hmask, h_nodal_in, rhs, uh_ice, vh_ice, dt)
  type(ice_shelf_dyn_CS), intent(in) :: CS
  type(ocean_grid_type),  intent(in) :: G
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in)        :: hmask
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(in)    :: h_nodal_in
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(out)   :: rhs
  real, dimension(SZDIB_(G),SZDJ_(G)),    intent(inout) :: uh_ice
  real, dimension(SZDI_(G),SZDJB_(G)),    intent(inout) :: vh_ice
  real,                                   intent(in)    :: dt !< RK-stage time step [T ~> s]
                                                              !! used to derive the per-face
                                                              !! CFL cap on the gated
                                                              !! artificial-viscosity coef.

  ! 2-point Gauss-Legendre on [0,1]
  real, parameter :: gp1 = 0.5 - 0.5/sqrt(3.0)
  real, parameter :: gp2 = 0.5 + 0.5/sqrt(3.0)
  real, parameter :: gw  = 0.5

  integer :: i, j, isc, iec, jsc, jec, qx, qy, gp, a, b
  real :: xi_q, eta_q, a_qp, d_qp
  real :: h_qp, u_qp, v_qp
  real :: N11, N21, N12, N22
  real :: dN_dxi_11, dN_dxi_21, dN_dxi_12, dN_dxi_22
  real :: dN_deta_11, dN_deta_21, dN_deta_12, dN_deta_22
  real :: dxCv_S, dxCv_N, dyCu_W, dyCu_E
  real :: t_face, t_co
  real :: u_at_qp, v_at_qp, h_upwind, flux_qp, face_flux_total
  real :: h_A_qp, h_B_qp     ! Face-QP corner thickness on the two sides of a DG face [Z ~> m]
  real :: u_mag_qp           ! Velocity magnitude sqrt(u^2+v^2) at a face QP [L T-1 ~> m s-1]
  real :: visc_flux_qp       ! Artificial-viscosity face flux per QP, antisymmetric [Z L2 T-1]
  real :: Hbar_A, Hbar_B     ! Cell-mean thickness on the two sides of a DG face [Z ~> m]
  real :: H_ref              ! Reference thickness for the smoothness ratio [Z ~> m]
  real :: r_face             ! Smoothness ratio |Delta h_eq|/H_ref [nondim]
  real :: sigma_face         ! Piecewise-linear smoothness ramp sigma(r_face) in [0,1] [nondim]
  real :: dx_perp            ! Across-face length scale at the current face [L ~> m]
  real :: coef_face          ! Per-face viscosity coefficient post smoothness gate [nondim]
  real :: u_eff_face_max     ! Max u_eff over the 2 face QPs, for the stability budget [L T-1 ~> m s-1]
  real :: amp_qp             ! Per-QP jump-mode rate amplification [nondim]
  real :: amp_face_max       ! Max amp_qp over the 2 face QPs [nondim]
  real :: ds_qp              ! Surface-elevation jump s_B - s_A at a face QP [Z ~> m]
  real :: ds_bar             ! Nonnegative mean-supported surface-jump allowance from
                             ! the cell means [Z ~> m]
  real :: ds_use             ! Surface jump driving the flux: full or excess [Z ~> m]
  real :: ds_face_max        ! Max |ds_use| over the 2 face QPs [Z ~> m]
  real :: ds_bar_face_max    ! Max |ds_bar| over the 2 face QPs (allowance diagnostic) [Z ~> m]
  real :: excess_frac_face   ! Max |ds_use|/|ds_qp| over the 2 face QPs (diagnostic) [nondim]
  real :: rate_face          ! Per-face semi-discrete diffusion rate c*u_eff_avg*ell/dx_perp [T-1]
  real :: scale_AB           ! Min of the two adjacent cells' per-cell CFL scale [nondim]
  real :: dh_eq              ! Well-balanced equivalent thickness jump (surface-jump-derived) [Z ~> m]
                             ! driving the viscosity flux.
  real :: bed_qp             ! Bed elevation at a face QP [Z ~> m]
  real :: rhoi_rhow_wb       ! Ice/ocean density ratio for the well-balanced jump [nondim]
  real :: u_floor_qp         ! Strain-rate-scaled velocity-independent floor at a face QP [L T-1 ~> m s-1]
  real :: u_eff_qp           ! |u_face| + u_floor for the viscosity flux and CFL cap [L T-1 ~> m s-1]
  real :: eps_e_face         ! Effective strain rate at the face midpoint [T-1]
  real :: dudx_f, dudy_f, dvdx_f, dvdy_f ! Face-midpoint velocity gradients [T-1]
  real :: u_mn, u_pl, v_mn, v_pl ! 4-corner-averaged cell velocities on the [L T-1 ~> m s-1]
                                  ! minus/plus side of the face
  real :: dh_eq_face_max     ! Max |Delta h_eq| over the 2 face QPs [Z ~> m]
  real :: u_mag_face_max     ! Max |u_mag_qp| over the 2 face QPs [L T-1 ~> m s-1]
  logical :: valid_A_visc, valid_B_visc ! Side A/B participates in DG art-visc face
                                         ! flux (hmask==1 interior, or hmask==3
                                         ! Dirichlet thickness BC used one-sided).
  ! The smoothness-gate thresholds (CS%dg_art_visc_r_lo/r_hi), baseline coefficient
  ! (CS%dg_art_visc_c_min), and per-cell jump-mode stability budget
  ! (CS%dg_art_visc_kcell) are runtime parameters; see their get_param descriptions.
  ! Pass-1 workspace for the two-pass per-cell CFL restructure of the viscosity branch:
  ! per-face desired coefficient, per-QP u_eff and Delta h_eq, and per-cell sum of face
  ! rates. Sized over the same index ranges as the surrounding face loops.
  real, dimension(SZDIB_(G),SZDJ_(G))   :: cK_E, rate_E
  real, dimension(SZDI_(G),SZDJB_(G))   :: cK_N, rate_N
  real, dimension(SZDIB_(G),SZDJ_(G),2) :: ueff_E, dheq_E
  real, dimension(SZDI_(G),SZDJB_(G),2) :: ueff_N, dheq_N
  logical, dimension(SZDIB_(G),SZDJ_(G)) :: active_E
  logical, dimension(SZDI_(G),SZDJB_(G)) :: active_N
  real, dimension(SZDI_(G),SZDJ_(G))    :: cell_scale
  real :: S_K                ! Per-cell sum of face rates for the CFL bound [T-1]
  integer :: i_lo, i_hi, j_lo, j_hi  ! One-sided clipping bounds for boundary strain
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: rhs_vol, rhs_face

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  rhoi_rhow_wb = CS%density_ice / CS%density_ocean_avg

  rhs(:,:,:,:) = 0.0
  rhs_vol(:,:,:,:) = 0.0
  rhs_face(:,:,:,:) = 0.0
  if (associated(CS%dg_art_visc_coef_u)) CS%dg_art_visc_coef_u(:,:) = 0.0
  if (associated(CS%dg_art_visc_coef_v)) CS%dg_art_visc_coef_v(:,:) = 0.0
  if (associated(CS%dg_art_visc_nu_u)) CS%dg_art_visc_nu_u(:,:) = 0.0
  if (associated(CS%dg_art_visc_nu_v)) CS%dg_art_visc_nu_v(:,:) = 0.0
  if (associated(CS%dg_art_visc_excess_frac_u)) CS%dg_art_visc_excess_frac_u(:,:) = 0.0
  if (associated(CS%dg_art_visc_excess_frac_v)) CS%dg_art_visc_excess_frac_v(:,:) = 0.0
  if (associated(CS%dg_art_visc_allow_u)) CS%dg_art_visc_allow_u(:,:) = 0.0
  if (associated(CS%dg_art_visc_allow_v)) CS%dg_art_visc_allow_v(:,:) = 0.0
  if (associated(CS%dg_slow_idle_face_u)) CS%dg_slow_idle_face_u(:,:) = 0.0
  if (associated(CS%dg_slow_idle_face_v)) CS%dg_slow_idle_face_v(:,:) = 0.0
  ! Reset to 1 (= unthrottled), not 0: a 0 would read as "fully throttled" in cells
  ! the cap loop never visits (non-ice cells, or the whole domain when the
  ! viscosity is disabled).
  if (associated(CS%dg_art_visc_cell_scale)) CS%dg_art_visc_cell_scale(:,:) = 1.0

  ! Volume integral over each cell.
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    dxCv_S = G%dxCv(i,j-1) ; dxCv_N = G%dxCv(i,j)
    dyCu_W = G%dyCu(i-1,j) ; dyCu_E = G%dyCu(i,j)

    do qy = 1, 2 ; do qx = 1, 2
      if (qx == 1) then ; xi_q  = gp1 ; else ; xi_q  = gp2 ; endif
      if (qy == 1) then ; eta_q = gp1 ; else ; eta_q = gp2 ; endif
      a_qp = dxCv_S*(1.0 - eta_q) + dxCv_N*eta_q
      d_qp = dyCu_W*(1.0 - xi_q)  + dyCu_E*xi_q

      N11 = (1.0-xi_q)*(1.0-eta_q) ; N21 = xi_q*(1.0-eta_q)
      N12 = (1.0-xi_q)*eta_q       ; N22 = xi_q*eta_q
      dN_dxi_11  = -(1.0 - eta_q) ; dN_dxi_21  =  (1.0 - eta_q)
      dN_dxi_12  = -eta_q         ; dN_dxi_22  =  eta_q
      dN_deta_11 = -(1.0 - xi_q)  ; dN_deta_21 = -xi_q
      dN_deta_12 =  (1.0 - xi_q)  ; dN_deta_22 =  xi_q

      h_qp = ((N11*h_nodal_in(i,j,1,1) + N22*h_nodal_in(i,j,2,2)) + &
              (N21*h_nodal_in(i,j,2,1) + N12*h_nodal_in(i,j,1,2)))
      u_qp = ((N11*CS%u_shelf(i-1,j-1) + N22*CS%u_shelf(i,j)) + &
              (N21*CS%u_shelf(i,j-1)   + N12*CS%u_shelf(i-1,j)))
      v_qp = ((N11*CS%v_shelf(i-1,j-1) + N22*CS%v_shelf(i,j)) + &
              (N21*CS%v_shelf(i,j-1)   + N12*CS%v_shelf(i-1,j)))

      ! Volume contribution at this QP for each test function N(a,b):
      ! + weight * h * ( u * dN/dxi * d + v * dN/deta * a )
      rhs_vol(i,j,1,1) = rhs_vol(i,j,1,1) + &
        gw*gw * h_qp * (u_qp * dN_dxi_11 * d_qp + v_qp * dN_deta_11 * a_qp)
      rhs_vol(i,j,2,1) = rhs_vol(i,j,2,1) + &
        gw*gw * h_qp * (u_qp * dN_dxi_21 * d_qp + v_qp * dN_deta_21 * a_qp)
      rhs_vol(i,j,1,2) = rhs_vol(i,j,1,2) + &
        gw*gw * h_qp * (u_qp * dN_dxi_12 * d_qp + v_qp * dN_deta_12 * a_qp)
      rhs_vol(i,j,2,2) = rhs_vol(i,j,2,2) + &
        gw*gw * h_qp * (u_qp * dN_dxi_22 * d_qp + v_qp * dN_deta_22 * a_qp)
    enddo ; enddo
  enddo ; enddo

  ! East-face fluxes between cells (i,j) and (i+1,j).
  do j = jsc, jec ; do i = isc-1, iec
    if (CS%u_face_mask(i,j) == 4.0) then
      face_flux_total = G%dyCu(i,j) * CS%u_flux_bdry_val(i,j)
      uh_ice(i,j) = uh_ice(i,j) + face_flux_total
      if (i >= isc .and. hmask(i,j) == 1.0) then
        rhs_face(i,j,2,1) = rhs_face(i,j,2,1) - 0.5*face_flux_total
        rhs_face(i,j,2,2) = rhs_face(i,j,2,2) - 0.5*face_flux_total
      endif
      if (i+1 <= iec .and. hmask(i+1,j) == 1.0) then
        rhs_face(i+1,j,1,1) = rhs_face(i+1,j,1,1) + 0.5*face_flux_total
        rhs_face(i+1,j,1,2) = rhs_face(i+1,j,1,2) + 0.5*face_flux_total
      endif
    else if (((i >= isc .and. (hmask(i,j) == 1.0 .or. hmask(i,j) == 3.0))) .or. &
             ((i+1 <= iec .and. (hmask(i+1,j) == 1.0 .or. hmask(i+1,j) == 3.0)))) then
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        u_at_qp = t_co*CS%u_shelf(i,j-1) + t_face*CS%u_shelf(i,j)
        if (u_at_qp >= 0.0) then
          if (hmask(i,j) == 3.0) then
            h_upwind = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
          elseif (hmask(i,j) == 1.0) then
            h_upwind = t_co*h_nodal_in(i,j,2,1) + t_face*h_nodal_in(i,j,2,2)
          else
            h_upwind = 0.0
          endif
        else
          if (hmask(i+1,j) == 3.0) then
            h_upwind = max(CS%h_bdry_val(i+1,j), CS%min_h_shelf)
          elseif (hmask(i+1,j) == 1.0) then
            h_upwind = t_co*h_nodal_in(i+1,j,1,1) + t_face*h_nodal_in(i+1,j,1,2)
          else
            h_upwind = 0.0
          endif
        endif
        h_upwind = max(h_upwind, 0.0)
        flux_qp = gw * u_at_qp * h_upwind * G%dyCu(i,j)
        uh_ice(i,j) = uh_ice(i,j) + flux_qp
        if (i >= isc .and. hmask(i,j) == 1.0) then
          rhs_face(i,j,2,1) = rhs_face(i,j,2,1) - flux_qp * t_co
          rhs_face(i,j,2,2) = rhs_face(i,j,2,2) - flux_qp * t_face
        endif
        if (i+1 <= iec .and. hmask(i+1,j) == 1.0) then
          rhs_face(i+1,j,1,1) = rhs_face(i+1,j,1,1) + flux_qp * t_co
          rhs_face(i+1,j,1,2) = rhs_face(i+1,j,1,2) + flux_qp * t_face
        endif
      enddo
    endif
  enddo ; enddo

  ! DG(1) artificial viscosity, face formulation.
  ! Per face, the antisymmetric flux is +gw * c_face * u_eff * Delta h_eq * ell_face
  ! at each QP, distributed conservatively to the two adjacent cells (-A, +B).
  ! Driver is the well-balanced equivalent thickness jump Delta h_eq derived from
  ! the surface jump, so a hydrostatically-continuous grounding line is not damped.
  ! With DG1_ART_VISC_EXCESS_JUMP, only the part of the surface jump in excess of
  ! the jump supported by the two cell-mean surfaces drives the flux, so standing
  ! mean-supported contrasts (e.g. shear margins) are not damped.
  ! c_face = c_min + (c_max - c_min) * sigma(r_face) ramps from c_min in smooth
  ! regions to c_max at shocks, where the gate ratio r_face is |Delta h_eq|/H_ref
  ! (legacy) or |[s]|/H_ref (DG1_ART_VISC_GATE_SURFACE): O(dx^2) smooth, O(1) at
  ! jumps. u_eff = advect_coef*|u_face| + strain_coef*eps_e*dx_perp covers both
  ! advective and deformation-driven excitation of the broken-Q1 mode. A per-cell
  ! SSP-RK2 stability budget (DG1_ART_VISC_KCELL/dt) on the summed jump-mode decay
  ! rates lambda_F = 4*amp*c*u_eff_max/dx_perp scales all faces of any cell that
  ! would exceed it; the factor 4*amp (>= 8) accounts for node localization (x2),
  ! the consistent-mass inverse at the face node (x2), and the flotation-branch
  ! flux Jacobian (amp = 2 uniform/harmonic, up to ~5.7 mixed-arithmetic) relative
  ! to the cell-mean rate. The excess map is 1-Lipschitz in the face traces, so
  ! this rate remains a valid bound with DG1_ART_VISC_EXCESS_JUMP. Weak one-sided
  ! imposition at Dirichlet thickness BCs (hmask==3): the existing hmask==1 write
  ! guards on rhs_face discard the BC-side update; the BC side uses h_bdry_val for
  ! h and Hbar.
  if (CS%dg_art_visc_c_max > 0.0) then
    active_E(:,:) = .false.
    cK_E(:,:)     = 0.0
    rate_E(:,:)   = 0.0
    ueff_E(:,:,:) = 0.0
    dheq_E(:,:,:) = 0.0
    active_N(:,:) = .false.
    cK_N(:,:)     = 0.0
    rate_N(:,:)   = 0.0
    ueff_N(:,:,:) = 0.0
    dheq_N(:,:,:) = 0.0
    cell_scale(:,:) = 1.0

    ! Pass 1, east faces: compute per-QP (u_eff, Delta h_eq), face-level c_face and
    ! semi-discrete diffusion rate, store for the scaled application in pass 2.
    do j = jsc, jec ; do i = isc-1, iec
      if (CS%u_face_mask(i,j) == 4.0) cycle  ! specified-flux face: viscosity undefined.
      valid_A_visc = (hmask(i,  j) == 1.0 .or. hmask(i,  j) == 3.0)
      valid_B_visc = (hmask(i+1,j) == 1.0 .or. hmask(i+1,j) == 3.0)
      if (.not. (valid_A_visc .and. valid_B_visc)) cycle
      if (.not. (hmask(i,j) == 1.0 .or. hmask(i+1,j) == 1.0)) cycle
      active_E(i,j) = .true.

      if (hmask(i,j) == 3.0) then
        Hbar_A = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
      else
        Hbar_A = nodal_cell_mean(h_nodal_in(i,  j,:,:), CS%cell_mean_w(i,  j,:,:))
      endif
      if (hmask(i+1,j) == 3.0) then
        Hbar_B = max(CS%h_bdry_val(i+1,j), CS%min_h_shelf)
      else
        Hbar_B = nodal_cell_mean(h_nodal_in(i+1,j,:,:), CS%cell_mean_w(i+1,j,:,:))
      endif
      H_ref = max(CS%min_h_shelf, 0.5*(Hbar_A + Hbar_B))
      dx_perp = G%dxCu(i,j)

      ! Boundary-aware one-sided cell-centred velocities for the across-face strain.
      ! When the i-1 or i+1 column is outside the data halo, collapse to the central
      ! column (first-order one-sided difference instead of zeroing the floor).
      i_lo = max(i-1, G%isd) ; i_hi = min(i+1, G%ied)
      if (CS%dg_art_visc_strain_coef > 0.0) then
        if (i_lo < i) then
          u_mn = 0.25*((CS%u_shelf(i_lo,j-1) + CS%u_shelf(i,j-1)) + &
                       (CS%u_shelf(i_lo,j  ) + CS%u_shelf(i,j  )))
          v_mn = 0.25*((CS%v_shelf(i_lo,j-1) + CS%v_shelf(i,j-1)) + &
                       (CS%v_shelf(i_lo,j  ) + CS%v_shelf(i,j  )))
        else
          u_mn = 0.5*(CS%u_shelf(i,j-1) + CS%u_shelf(i,j))
          v_mn = 0.5*(CS%v_shelf(i,j-1) + CS%v_shelf(i,j))
        endif
        if (i_hi > i) then
          u_pl = 0.25*((CS%u_shelf(i   ,j-1) + CS%u_shelf(i_hi,j-1)) + &
                       (CS%u_shelf(i   ,j  ) + CS%u_shelf(i_hi,j  )))
          v_pl = 0.25*((CS%v_shelf(i   ,j-1) + CS%v_shelf(i_hi,j-1)) + &
                       (CS%v_shelf(i   ,j  ) + CS%v_shelf(i_hi,j  )))
        else
          u_pl = 0.5*(CS%u_shelf(i,j-1) + CS%u_shelf(i,j))
          v_pl = 0.5*(CS%v_shelf(i,j-1) + CS%v_shelf(i,j))
        endif
        dudx_f = (u_pl - u_mn) / dx_perp
        dvdx_f = (v_pl - v_mn) / dx_perp
        dudy_f = (CS%u_shelf(i,j) - CS%u_shelf(i,j-1)) / G%dyCu(i,j)
        dvdy_f = (CS%v_shelf(i,j) - CS%v_shelf(i,j-1)) / G%dyCu(i,j)
        eps_e_face = dg1_face_eps_eff(dudx_f, dudy_f, dvdx_f, dvdy_f)
      else
        eps_e_face = 0.0
      endif

      dh_eq_face_max = 0.0
      u_mag_face_max = 0.0
      u_eff_face_max = 0.0
      amp_face_max = 0.0
      ds_face_max = 0.0
      ds_bar_face_max = 0.0
      excess_frac_face = 0.0
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        u_at_qp = t_co*CS%u_shelf(i,j-1) + t_face*CS%u_shelf(i,j)
        v_at_qp = t_co*CS%v_shelf(i,j-1) + t_face*CS%v_shelf(i,j)
        u_mag_qp = sqrt(u_at_qp*u_at_qp + v_at_qp*v_at_qp)
        if (hmask(i,j) == 3.0) then
          h_A_qp = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
        else
          h_A_qp = t_co*h_nodal_in(i,  j,2,1) + t_face*h_nodal_in(i,  j,2,2)
        endif
        if (hmask(i+1,j) == 3.0) then
          h_B_qp = max(CS%h_bdry_val(i+1,j), CS%min_h_shelf)
        else
          h_B_qp = t_co*h_nodal_in(i+1,j,1,1) + t_face*h_nodal_in(i+1,j,1,2)
        endif
        bed_qp = t_co*CS%bed_node(i,j-1) + t_face*CS%bed_node(i,j)
        if (CS%dg_art_visc_excess_jump) then
          ds_qp  = dg1_wb_surface_jump(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb)
          ! Mean-supported allowance: per-side surfaces evaluated from the cell
          ! means at the face-QP bed; with BRANCH_MAX, the max over the admissible
          ! flotation-branch assignments (see dg1_wb_mean_allowance).
          ds_bar = dg1_wb_mean_allowance(Hbar_A, Hbar_B, h_A_qp, h_B_qp, bed_qp, &
                                         rhoi_rhow_wb, CS%dg_art_visc_excess_branch_max)
          ! Excess over the supported band [-ds_bar, +ds_bar].
          ds_use = ds_qp - min(max(ds_qp, -ds_bar), ds_bar)
          dh_eq = ds_use * dg1_wb_slope_mean(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb, &
                                             CS%dg_art_visc_wb_harmonic)
          ds_bar_face_max = max(ds_bar_face_max, ds_bar)
          if (abs(ds_qp) > 0.0) &
            excess_frac_face = max(excess_frac_face, abs(ds_use) / abs(ds_qp))
        else
          if (CS%dg_art_visc_gate_surface) &
            ds_use = dg1_wb_surface_jump(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb)
          dh_eq = dg1_wb_equiv_jump(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb, &
                                    CS%dg_art_visc_wb_harmonic)
        endif
        amp_qp = dg1_wb_jump_rate_amp(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb, &
                                      CS%dg_art_visc_wb_harmonic)
        u_floor_qp = CS%dg_art_visc_strain_coef * eps_e_face * dx_perp
        u_eff_qp = CS%dg_art_visc_advect_coef * u_mag_qp + u_floor_qp

        ueff_E(i,j,gp) = u_eff_qp
        dheq_E(i,j,gp) = dh_eq

        dh_eq_face_max = max(dh_eq_face_max, abs(dh_eq))
        if (CS%dg_art_visc_gate_surface) ds_face_max = max(ds_face_max, abs(ds_use))
        amp_face_max = max(amp_face_max, amp_qp)
        u_eff_face_max = max(u_eff_face_max, u_eff_qp)
        u_mag_face_max = max(u_mag_face_max, u_mag_qp)
      enddo

      if (CS%dg_art_visc_excess_jump) then
        if (associated(CS%dg_art_visc_excess_frac_u)) &
          CS%dg_art_visc_excess_frac_u(i,j) = excess_frac_face
        if (associated(CS%dg_art_visc_allow_u)) &
          CS%dg_art_visc_allow_u(i,j) = ds_bar_face_max
      endif

      if (CS%dg_art_visc_gate_surface) then
        ! Surface-cliff gate: R_LO/R_HI are the cliff heights, as fractions of the
        ! local mean thickness, at which damping starts and saturates.
        r_face = ds_face_max / H_ref
      else
        r_face = dh_eq_face_max / H_ref
      endif
      if (r_face <= CS%dg_art_visc_r_lo) then
        sigma_face = 0.0
      elseif (r_face >= CS%dg_art_visc_r_hi) then
        sigma_face = 1.0
      else
        sigma_face = (r_face - CS%dg_art_visc_r_lo) / &
                     (CS%dg_art_visc_r_hi - CS%dg_art_visc_r_lo)
      endif
      coef_face = CS%dg_art_visc_c_min + &
                  (CS%dg_art_visc_c_max - CS%dg_art_visc_c_min) * sigma_face
      cK_E(i,j) = coef_face
      ! Per-face semi-discrete decay rate of the face-node jump mode,
      ! lambda_F = 4 * amp * c * u_eff / dx_perp. The factor 4*amp (>= 8) accounts
      ! for node localization (x2), the consistent-mass inverse at the face node
      ! (x2), and the flotation-branch flux Jacobian (amp = 2 uniform/harmonic, up
      ! to ~5.7 mixed-arithmetic), relative to the cell-mean rate c*u_eff/dx_perp.
      ! Max-over-QP u_eff bounds the worst QP. No ell factor (it cancels between
      ! the face-integrated flux and the cell area).
      rate_face = 4.0 * amp_face_max * coef_face * u_eff_face_max / dx_perp
      rate_E(i,j) = rate_face

      ! Stagnant-jump diagnostic: |u| and eps_e both tiny while the equivalent jump
      ! is non-negligible. Mark with 1.0 (otherwise stays at 0 from the reset).
      if (associated(CS%dg_slow_idle_face_u)) then
        if (u_mag_face_max < CS%dg_slow_idle_u_tiny .and. &
            eps_e_face     < CS%dg_slow_idle_eps_tiny .and. &
            dh_eq_face_max > CS%dg_slow_idle_s_tol) then
          CS%dg_slow_idle_face_u(i,j) = 1.0
        endif
      endif
    enddo ; enddo
  endif

  ! North-face fluxes between cells (i,j) and (i,j+1).
  do j = jsc-1, jec ; do i = isc, iec
    if (CS%v_face_mask(i,j) == 4.0) then
      face_flux_total = G%dxCv(i,j) * CS%v_flux_bdry_val(i,j)
      vh_ice(i,j) = vh_ice(i,j) + face_flux_total
      if (j >= jsc .and. hmask(i,j) == 1.0) then
        rhs_face(i,j,1,2) = rhs_face(i,j,1,2) - 0.5*face_flux_total
        rhs_face(i,j,2,2) = rhs_face(i,j,2,2) - 0.5*face_flux_total
      endif
      if (j+1 <= jec .and. hmask(i,j+1) == 1.0) then
        rhs_face(i,j+1,1,1) = rhs_face(i,j+1,1,1) + 0.5*face_flux_total
        rhs_face(i,j+1,2,1) = rhs_face(i,j+1,2,1) + 0.5*face_flux_total
      endif
    else if (((j >= jsc .and. (hmask(i,j) == 1.0 .or. hmask(i,j) == 3.0))) .or. &
             ((j+1 <= jec .and. (hmask(i,j+1) == 1.0 .or. hmask(i,j+1) == 3.0)))) then
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        v_at_qp = t_co*CS%v_shelf(i-1,j) + t_face*CS%v_shelf(i,j)
        if (v_at_qp >= 0.0) then
          if (hmask(i,j) == 3.0) then
            h_upwind = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
          elseif (hmask(i,j) == 1.0) then
            h_upwind = t_co*h_nodal_in(i,j,1,2) + t_face*h_nodal_in(i,j,2,2)
          else
            h_upwind = 0.0
          endif
        else
          if (hmask(i,j+1) == 3.0) then
            h_upwind = max(CS%h_bdry_val(i,j+1), CS%min_h_shelf)
          elseif (hmask(i,j+1) == 1.0) then
            h_upwind = t_co*h_nodal_in(i,j+1,1,1) + t_face*h_nodal_in(i,j+1,2,1)
          else
            h_upwind = 0.0
          endif
        endif
        h_upwind = max(h_upwind, 0.0)
        flux_qp = gw * v_at_qp * h_upwind * G%dxCv(i,j)
        vh_ice(i,j) = vh_ice(i,j) + flux_qp
        if (j >= jsc .and. hmask(i,j) == 1.0) then
          rhs_face(i,j,1,2) = rhs_face(i,j,1,2) - flux_qp * t_co
          rhs_face(i,j,2,2) = rhs_face(i,j,2,2) - flux_qp * t_face
        endif
        if (j+1 <= jec .and. hmask(i,j+1) == 1.0) then
          rhs_face(i,j+1,1,1) = rhs_face(i,j+1,1,1) + flux_qp * t_co
          rhs_face(i,j+1,2,1) = rhs_face(i,j+1,2,1) + flux_qp * t_face
        endif
      enddo
    endif
  enddo ; enddo

  ! Pass 1, north faces (analogous to east-face block above).
  if (CS%dg_art_visc_c_max > 0.0) then
    do j = jsc-1, jec ; do i = isc, iec
      if (CS%v_face_mask(i,j) == 4.0) cycle
      valid_A_visc = (hmask(i,j  ) == 1.0 .or. hmask(i,j  ) == 3.0)
      valid_B_visc = (hmask(i,j+1) == 1.0 .or. hmask(i,j+1) == 3.0)
      if (.not. (valid_A_visc .and. valid_B_visc)) cycle
      if (.not. (hmask(i,j) == 1.0 .or. hmask(i,j+1) == 1.0)) cycle
      active_N(i,j) = .true.

      if (hmask(i,j) == 3.0) then
        Hbar_A = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
      else
        Hbar_A = nodal_cell_mean(h_nodal_in(i,j,  :,:), CS%cell_mean_w(i,j,  :,:))
      endif
      if (hmask(i,j+1) == 3.0) then
        Hbar_B = max(CS%h_bdry_val(i,j+1), CS%min_h_shelf)
      else
        Hbar_B = nodal_cell_mean(h_nodal_in(i,j+1,:,:), CS%cell_mean_w(i,j+1,:,:))
      endif
      H_ref = max(CS%min_h_shelf, 0.5*(Hbar_A + Hbar_B))
      dx_perp = G%dyCv(i,j)

      j_lo = max(j-1, G%jsd) ; j_hi = min(j+1, G%jed)
      if (CS%dg_art_visc_strain_coef > 0.0) then
        if (j_lo < j) then
          u_mn = 0.25*((CS%u_shelf(i-1,j_lo) + CS%u_shelf(i,j_lo)) + &
                       (CS%u_shelf(i-1,j   ) + CS%u_shelf(i,j   )))
          v_mn = 0.25*((CS%v_shelf(i-1,j_lo) + CS%v_shelf(i,j_lo)) + &
                       (CS%v_shelf(i-1,j   ) + CS%v_shelf(i,j   )))
        else
          u_mn = 0.5*(CS%u_shelf(i-1,j) + CS%u_shelf(i,j))
          v_mn = 0.5*(CS%v_shelf(i-1,j) + CS%v_shelf(i,j))
        endif
        if (j_hi > j) then
          u_pl = 0.25*((CS%u_shelf(i-1,j   ) + CS%u_shelf(i,j   )) + &
                       (CS%u_shelf(i-1,j_hi) + CS%u_shelf(i,j_hi)))
          v_pl = 0.25*((CS%v_shelf(i-1,j   ) + CS%v_shelf(i,j   )) + &
                       (CS%v_shelf(i-1,j_hi) + CS%v_shelf(i,j_hi)))
        else
          u_pl = 0.5*(CS%u_shelf(i-1,j) + CS%u_shelf(i,j))
          v_pl = 0.5*(CS%v_shelf(i-1,j) + CS%v_shelf(i,j))
        endif
        dudy_f = (u_pl - u_mn) / dx_perp
        dvdy_f = (v_pl - v_mn) / dx_perp
        dudx_f = (CS%u_shelf(i,j) - CS%u_shelf(i-1,j)) / G%dxCv(i,j)
        dvdx_f = (CS%v_shelf(i,j) - CS%v_shelf(i-1,j)) / G%dxCv(i,j)
        eps_e_face = dg1_face_eps_eff(dudx_f, dudy_f, dvdx_f, dvdy_f)
      else
        eps_e_face = 0.0
      endif

      dh_eq_face_max = 0.0
      u_mag_face_max = 0.0
      u_eff_face_max = 0.0
      amp_face_max = 0.0
      ds_face_max = 0.0
      ds_bar_face_max = 0.0
      excess_frac_face = 0.0
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        u_at_qp = t_co*CS%u_shelf(i-1,j) + t_face*CS%u_shelf(i,j)
        v_at_qp = t_co*CS%v_shelf(i-1,j) + t_face*CS%v_shelf(i,j)
        u_mag_qp = sqrt(u_at_qp*u_at_qp + v_at_qp*v_at_qp)
        if (hmask(i,j) == 3.0) then
          h_A_qp = max(CS%h_bdry_val(i,j), CS%min_h_shelf)
        else
          h_A_qp = t_co*h_nodal_in(i,j,  1,2) + t_face*h_nodal_in(i,j,  2,2)
        endif
        if (hmask(i,j+1) == 3.0) then
          h_B_qp = max(CS%h_bdry_val(i,j+1), CS%min_h_shelf)
        else
          h_B_qp = t_co*h_nodal_in(i,j+1,1,1) + t_face*h_nodal_in(i,j+1,2,1)
        endif
        bed_qp = t_co*CS%bed_node(i-1,j) + t_face*CS%bed_node(i,j)
        if (CS%dg_art_visc_excess_jump) then
          ds_qp  = dg1_wb_surface_jump(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb)
          ds_bar = dg1_wb_mean_allowance(Hbar_A, Hbar_B, h_A_qp, h_B_qp, bed_qp, &
                                         rhoi_rhow_wb, CS%dg_art_visc_excess_branch_max)
          ds_use = ds_qp - min(max(ds_qp, -ds_bar), ds_bar)
          dh_eq = ds_use * dg1_wb_slope_mean(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb, &
                                             CS%dg_art_visc_wb_harmonic)
          ds_bar_face_max = max(ds_bar_face_max, ds_bar)
          if (abs(ds_qp) > 0.0) &
            excess_frac_face = max(excess_frac_face, abs(ds_use) / abs(ds_qp))
        else
          if (CS%dg_art_visc_gate_surface) &
            ds_use = dg1_wb_surface_jump(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb)
          dh_eq = dg1_wb_equiv_jump(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb, &
                                    CS%dg_art_visc_wb_harmonic)
        endif
        amp_qp = dg1_wb_jump_rate_amp(h_A_qp, h_B_qp, bed_qp, rhoi_rhow_wb, &
                                      CS%dg_art_visc_wb_harmonic)
        u_floor_qp = CS%dg_art_visc_strain_coef * eps_e_face * dx_perp
        u_eff_qp = CS%dg_art_visc_advect_coef * u_mag_qp + u_floor_qp

        ueff_N(i,j,gp) = u_eff_qp
        dheq_N(i,j,gp) = dh_eq

        dh_eq_face_max = max(dh_eq_face_max, abs(dh_eq))
        if (CS%dg_art_visc_gate_surface) ds_face_max = max(ds_face_max, abs(ds_use))
        amp_face_max = max(amp_face_max, amp_qp)
        u_eff_face_max = max(u_eff_face_max, u_eff_qp)
        u_mag_face_max = max(u_mag_face_max, u_mag_qp)
      enddo

      if (CS%dg_art_visc_excess_jump) then
        if (associated(CS%dg_art_visc_excess_frac_v)) &
          CS%dg_art_visc_excess_frac_v(i,j) = excess_frac_face
        if (associated(CS%dg_art_visc_allow_v)) &
          CS%dg_art_visc_allow_v(i,j) = ds_bar_face_max
      endif

      if (CS%dg_art_visc_gate_surface) then
        r_face = ds_face_max / H_ref
      else
        r_face = dh_eq_face_max / H_ref
      endif
      if (r_face <= CS%dg_art_visc_r_lo) then
        sigma_face = 0.0
      elseif (r_face >= CS%dg_art_visc_r_hi) then
        sigma_face = 1.0
      else
        sigma_face = (r_face - CS%dg_art_visc_r_lo) / &
                     (CS%dg_art_visc_r_hi - CS%dg_art_visc_r_lo)
      endif
      coef_face = CS%dg_art_visc_c_min + &
                  (CS%dg_art_visc_c_max - CS%dg_art_visc_c_min) * sigma_face
      cK_N(i,j) = coef_face
      rate_face = 4.0 * amp_face_max * coef_face * u_eff_face_max / dx_perp
      rate_N(i,j) = rate_face

      if (associated(CS%dg_slow_idle_face_v)) then
        if (u_mag_face_max < CS%dg_slow_idle_u_tiny .and. &
            eps_e_face     < CS%dg_slow_idle_eps_tiny .and. &
            dh_eq_face_max > CS%dg_slow_idle_s_tol) then
          CS%dg_slow_idle_face_v(i,j) = 1.0
        endif
      endif
    enddo ; enddo

    ! Aggregate per-cell rate sum and form scale = min(1, kcell/(S_K*dt)). Each PE
    ! computes scales for its owned cells (all four face rates of an owned cell are
    ! available locally); non-ice and BC cells keep scale 1 so the min(scale_A,
    ! scale_B) below picks up the interior cell's budget at one-sided faces.
    do j = jsc, jec ; do i = isc, iec
      if (hmask(i,j) /= 1.0) cycle
      S_K = 0.0
      if (active_E(i-1,j)) S_K = S_K + rate_E(i-1,j)
      if (active_E(i  ,j)) S_K = S_K + rate_E(i  ,j)
      if (active_N(i,j-1)) S_K = S_K + rate_N(i,j-1)
      if (active_N(i,j  )) S_K = S_K + rate_N(i,j  )
      if (S_K*dt > CS%dg_art_visc_kcell) then
        cell_scale(i,j) = CS%dg_art_visc_kcell / (S_K*dt)
      else
        cell_scale(i,j) = 1.0
      endif
      if (associated(CS%dg_art_visc_cell_scale)) &
        CS%dg_art_visc_cell_scale(i,j) = cell_scale(i,j)
    enddo ; enddo
    ! Halo-update the scales so that a face on a PE boundary (or a reentrant seam)
    ! sees both adjacent cells' budgets: without this, each side would apply only
    ! its own cell's throttle, breaking the antisymmetric flux pair (local
    ! non-conservation) and making an engaged cap layout-dependent.
    call pass_var(cell_scale, G%domain)

    ! Pass 2, east faces: scaled application. Per-face scale_AB is the min of the
    ! two adjacent cells' scales; after the halo update both sides are valid on
    ! every active face (non-ice and BC cells carry scale 1), so the same value is
    ! formed on whichever PE computes the face and the antisymmetric pair is exact.
    do j = jsc, jec ; do i = isc-1, iec
      if (.not. active_E(i,j)) cycle
      scale_AB = min(cell_scale(i,j), cell_scale(i+1,j))
      coef_face = cK_E(i,j) * scale_AB
      if (associated(CS%dg_art_visc_coef_u)) CS%dg_art_visc_coef_u(i,j) = coef_face
      if (associated(CS%dg_art_visc_nu_u)) &
        CS%dg_art_visc_nu_u(i,j) = coef_face * 0.5*(ueff_E(i,j,1) + ueff_E(i,j,2)) * G%dxCu(i,j)
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        visc_flux_qp = gw * coef_face * ueff_E(i,j,gp) * dheq_E(i,j,gp) * G%dyCu(i,j)
        if (i >= isc .and. hmask(i,j) == 1.0) then
          rhs_face(i,j,2,1) = rhs_face(i,j,2,1) + visc_flux_qp * t_co
          rhs_face(i,j,2,2) = rhs_face(i,j,2,2) + visc_flux_qp * t_face
        endif
        if (i+1 <= iec .and. hmask(i+1,j) == 1.0) then
          rhs_face(i+1,j,1,1) = rhs_face(i+1,j,1,1) - visc_flux_qp * t_co
          rhs_face(i+1,j,1,2) = rhs_face(i+1,j,1,2) - visc_flux_qp * t_face
        endif
      enddo
    enddo ; enddo

    ! Pass 2, north faces.
    do j = jsc-1, jec ; do i = isc, iec
      if (.not. active_N(i,j)) cycle
      scale_AB = min(cell_scale(i,j), cell_scale(i,j+1))
      coef_face = cK_N(i,j) * scale_AB
      if (associated(CS%dg_art_visc_coef_v)) CS%dg_art_visc_coef_v(i,j) = coef_face
      if (associated(CS%dg_art_visc_nu_v)) &
        CS%dg_art_visc_nu_v(i,j) = coef_face * 0.5*(ueff_N(i,j,1) + ueff_N(i,j,2)) * G%dyCv(i,j)
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        visc_flux_qp = gw * coef_face * ueff_N(i,j,gp) * dheq_N(i,j,gp) * G%dxCv(i,j)
        if (j >= jsc .and. hmask(i,j) == 1.0) then
          rhs_face(i,j,1,2) = rhs_face(i,j,1,2) + visc_flux_qp * t_co
          rhs_face(i,j,2,2) = rhs_face(i,j,2,2) + visc_flux_qp * t_face
        endif
        if (j+1 <= jec .and. hmask(i,j+1) == 1.0) then
          rhs_face(i,j+1,1,1) = rhs_face(i,j+1,1,1) - visc_flux_qp * t_co
          rhs_face(i,j+1,2,1) = rhs_face(i,j+1,2,1) - visc_flux_qp * t_face
        endif
      enddo
    enddo ; enddo
  endif

  ! Reduce volume + face.
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    do b = 1, 2 ; do a = 1, 2
      rhs(i,j,a,b) = rhs_vol(i,j,a,b) + rhs_face(i,j,a,b)
    enddo ; enddo
  enddo ; enddo
end subroutine DG1_nodal_spatial_operator

!> Advect h_nodal one time step with SSP-RK2 + Barth-Jespersen + Liu positivity.
!! Operates directly on CS%h_nodal (mutates the authoritative nodal storage).
subroutine ice_shelf_advect_DG1_nodal(CS, ISS, G, time_step, hmask, uh_ice, vh_ice)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ice_shelf_state),  intent(in)    :: ISS
  type(ocean_grid_type),  intent(inout) :: G
  real,                   intent(in)    :: time_step
  real, dimension(SZDI_(G),SZDJ_(G)),       intent(inout) :: hmask
  real, dimension(SZDIB_(G),SZDJ_(G)),      intent(inout) :: uh_ice
  real, dimension(SZDI_(G),SZDJB_(G)),      intent(inout) :: vh_ice

  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: h0, h_curr, rhs
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: S_node ! Q1 nodal source per cell [Z T-1].
                                                   ! Continuous across cell faces by
                                                   ! construction; integrating it against
                                                   ! the local mass matrix recovers the
                                                   ! source contribution to dh_nodal/dt.
  real, dimension(2,2) :: dh
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed, a, b

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  uh_ice(:,:) = 0.0 ; vh_ice(:,:) = 0.0

  ! Dirichlet thickness BC: snap all 4 corners of BC cells to h_bdry_val.
  do j = jsd, jed ; do i = isd, ied
    if (CS%h_bdry_val(i,j) /= 0.0) CS%h_nodal(i,j,:,:) = CS%h_bdry_val(i,j)
  enddo ; enddo
  h0(:,:,:,:) = CS%h_nodal(:,:,:,:)
  call pass_corner_field(h0, G)

  ! Project the cell-mean source rate accumulated since the last advect
  ! (basal melt + surface SMB, units Z T-1) onto a continuous Q1 nodal field.
  ! For a Q1-projected source, the consistent Galerkin treatment M*dh/dt = -L
  ! + M*S_node collapses after M^-1 to dh/dt += S_node element-wise, so the
  ! source enters each SSP-RK2 stage as a simple additive term on dh.
  call project_h_source_rate_to_nodes(CS, ISS, G, S_node)

  ! Stage 1: positivity floor -> hierarchical limiter -> spatial op -> M^-1 -> Euler step.
  if (CS%nodal_positivity) call nodal_positivity_limit(CS, G, ISS)
  if (CS%dg_hierarchical_lim) call nodal_surface_slope_limit(CS, G, ISS)
  call pass_corner_field(CS%h_nodal, G)
  call DG1_nodal_spatial_operator(CS, G, hmask, CS%h_nodal, rhs, uh_ice, vh_ice, time_step)
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    call apply_nodal_DG_mass_inverse(CS%Minv_xi(i,j,:,:), CS%Minv_eta(i,j,:,:), &
                                     rhs(i,j,:,:), dh)
    do b = 1, 2 ; do a = 1, 2
      CS%h_nodal(i,j,a,b) = h0(i,j,a,b) + time_step * (dh(a,b) + S_node(i,j,a,b))
    enddo ; enddo
  enddo ; enddo
  call pass_corner_field(CS%h_nodal, G)

  ! Stage 2: positivity floor -> hierarchical limiter -> spatial op -> M^-1 -> SSP-RK2 combine.
  if (CS%nodal_positivity) call nodal_positivity_limit(CS, G, ISS)
  if (CS%dg_hierarchical_lim) call nodal_surface_slope_limit(CS, G, ISS)
  h_curr(:,:,:,:) = CS%h_nodal(:,:,:,:)
  call DG1_nodal_spatial_operator(CS, G, hmask, h_curr, rhs, uh_ice, vh_ice, time_step)
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    call apply_nodal_DG_mass_inverse(CS%Minv_xi(i,j,:,:), CS%Minv_eta(i,j,:,:), &
                                     rhs(i,j,:,:), dh)
    do b = 1, 2 ; do a = 1, 2
      CS%h_nodal(i,j,a,b) = 0.5*h0(i,j,a,b) &
                           + 0.5*(h_curr(i,j,a,b) + time_step*(dh(a,b) + S_node(i,j,a,b)))
    enddo ; enddo
  enddo ; enddo
  call pass_corner_field(CS%h_nodal, G)

  ! Snapshot the rate just consumed for the h_source_rate diagnostic, then
  ! reset the buffer so subsequent melt/SMB callers start from a clean slate.
  CS%h_source_rate_last(:,:) = CS%h_source_rate(:,:)
  CS%h_source_rate(:,:) = 0.0

  ! Final positivity floor + optional hierarchical limit on the SSP-RK2 result.
  if (CS%nodal_positivity) call nodal_positivity_limit(CS, G, ISS)
  if (CS%dg_hierarchical_lim) call nodal_surface_slope_limit(CS, G, ISS)
  call pass_corner_field(CS%h_nodal, G)

  ! Average uh_ice, vh_ice over the 2 stages (SSP-RK2 equal weight).
  uh_ice(:,:) = 0.5 * uh_ice(:,:)
  vh_ice(:,:) = 0.5 * vh_ice(:,:)
  call pass_vector(uh_ice, vh_ice, G%domain, TO_ALL, CGRID_NE)
end subroutine ice_shelf_advect_DG1_nodal

end module MOM_ice_shelf_dynamics
