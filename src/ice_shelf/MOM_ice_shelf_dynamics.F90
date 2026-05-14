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
implicit none ; private

#include <MOM_memory.h>

public register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn, update_ice_shelf, IS_dynamics_post_data
public ice_time_step_CFL, ice_shelf_dyn_end, change_in_draft, write_ice_shelf_energy
public shelf_advance_front, ice_shelf_min_thickness_calve, calve_to_mask, volume_above_floatation
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
  real, pointer, dimension(:,:) :: h_x => NULL()       !< DG(1) x-slope moment of ice thickness [Z ~> m].
                                                       !! h(xi,eta) = h_shelf + h_x*xi + h_y*eta, where
                                                       !! (xi,eta) in [-0.5, 0.5] are local cell coordinates.
  real, pointer, dimension(:,:) :: h_y => NULL()       !< DG(1) y-slope moment of ice thickness [Z ~> m].

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
  integer :: dg1_limiter_choice   !< Slope-limiter choice for DG(1) thickness:
                                  !! 0 = none, 1 = minmod (legacy), 2 = Venkatakrishnan.
  real :: dg1_limiter_M           !< TVB-style curvature bound for the Venkatakrishnan
                                  !! smooth-extremum protection band: eps = M * dx_local^2
                                  !! [Z L-2 ~> m-1].
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
             id_bed_node = -1, id_h_x = -1, id_h_y = -1, &
             id_phi_x_DG = -1, id_phi_y_DG = -1, &
             id_phi_x_FV = -1, id_phi_y_FV = -1
  real, pointer, dimension(:,:) :: phi_x_DG => NULL() !< DG(1) limiter factor in x at last advection
                                                       !! call [nondim], in [0,1].
  real, pointer, dimension(:,:) :: phi_y_DG => NULL() !< DG(1) limiter factor in y at last advection
                                                       !! call [nondim], in [0,1].
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
    allocate(CS%h_x(isd:ied,jsd:jed), source=0.0)
    allocate(CS%h_y(isd:ied,jsd:jed), source=0.0)
    allocate(CS%phi_x_DG(isd:ied,jsd:jed), source=1.0)
    allocate(CS%phi_y_DG(isd:ied,jsd:jed), source=1.0)
    allocate(CS%phi_x_FV(IsdB:IedB,jsd:jed), source=1.0)
    allocate(CS%phi_y_FV(isd:ied,JsdB:JedB), source=1.0)
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
    call register_restart_field(CS%h_x, "h_x_DG", .true., restart_CS, &
                                "DG(1) x-slope moment of ice thickness", "m", conversion=US%Z_to_m)
    call register_restart_field(CS%h_y, "h_y_DG", .true., restart_CS, &
                                "DG(1) y-slope moment of ice thickness", "m", conversion=US%Z_to_m)
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
    call read_dg1_limiter_params(param_file, mdl, CS, US)
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
        ! DG(1) cold-start slope init: seed from neighbour cell-mean central
        ! differences, then apply the slope limiter. Seeding (rather than
        ! zeroing) is required because DG1_slope_limit uses minmod3 of the
        ! existing slope with neighbour differences and minmod3 returns zero
        ! whenever any argument is zero, so a zero seed would persist.
        ! On restart h_x,h_y come from the restart file and this branch is skipped.
        call pass_var(ISS%h_shelf, G%domain)
        CS%h_x(:,:) = 0.0 ; CS%h_y(:,:) = 0.0
        do j=G%jsc,G%jec ; do i=G%isc,G%iec
          if (ISS%hmask(i,j) == 1) then
            valid_E = (ISS%hmask(i+1,j) == 1 .or. ISS%hmask(i+1,j) == 3)
            valid_W = (ISS%hmask(i-1,j) == 1 .or. ISS%hmask(i-1,j) == 3)
            valid_N = (ISS%hmask(i,j+1) == 1 .or. ISS%hmask(i,j+1) == 3)
            valid_S = (ISS%hmask(i,j-1) == 1 .or. ISS%hmask(i,j-1) == 3)
            ! For hmask==3 neighbours, h_bdry_val is the FACE value at distance
            ! dx/2 (see ice_shelf_advect_thickness_x). Mirror through the face
            ! to get an effective cell-mean at distance dx for centred FD.
            if (ISS%hmask(i+1,j) == 3) then
              h_E = 2.0 * CS%h_bdry_val(i+1,j) - ISS%h_shelf(i,j)
            else
              h_E = ISS%h_shelf(i+1,j)
            endif
            if (ISS%hmask(i-1,j) == 3) then
              h_W = 2.0 * CS%h_bdry_val(i-1,j) - ISS%h_shelf(i,j)
            else
              h_W = ISS%h_shelf(i-1,j)
            endif
            if (ISS%hmask(i,j+1) == 3) then
              h_N = 2.0 * CS%h_bdry_val(i,j+1) - ISS%h_shelf(i,j)
            else
              h_N = ISS%h_shelf(i,j+1)
            endif
            if (ISS%hmask(i,j-1) == 3) then
              h_S = 2.0 * CS%h_bdry_val(i,j-1) - ISS%h_shelf(i,j)
            else
              h_S = ISS%h_shelf(i,j-1)
            endif
            ! Centred difference where both neighbours are valid; one-sided
            ! difference at the ice front so the slope is preserved.
            if (valid_E .and. valid_W) then
              CS%h_x(i,j) = 0.5 * (h_E - h_W)
            elseif (valid_W) then
              CS%h_x(i,j) = ISS%h_shelf(i,j) - h_W
            elseif (valid_E) then
              CS%h_x(i,j) = h_E - ISS%h_shelf(i,j)
            endif
            if (valid_N .and. valid_S) then
              CS%h_y(i,j) = 0.5 * (h_N - h_S)
            elseif (valid_S) then
              CS%h_y(i,j) = ISS%h_shelf(i,j) - h_S
            elseif (valid_N) then
              CS%h_y(i,j) = h_N - ISS%h_shelf(i,j)
            endif
          endif
        enddo ; enddo
        call DG1_slope_limit(G, ISS%h_shelf, CS%h_x, CS%h_y, ISS%hmask, CS%h_bdry_val, &
                             CS%dg1_limiter_choice, CS%dg1_limiter_M)
        call pass_vector(CS%h_x, CS%h_y, G%domain, TO_ALL, AGRID)
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
      CS%id_h_x = register_diag_field('ice_shelf_model','h_x',CS%diag%axesT1, Time, &
         'DG(1) x-slope moment of ice thickness (h = hbar + h_x*xi + h_y*eta)', &
         'm', conversion=US%Z_to_m)
      CS%id_h_y = register_diag_field('ice_shelf_model','h_y',CS%diag%axesT1, Time, &
         'DG(1) y-slope moment of ice thickness', 'm', conversion=US%Z_to_m)
      CS%id_phi_x_DG = register_diag_field('ice_shelf_model','phi_x_DG',CS%diag%axesT1, Time, &
         'DG(1) slope-limiter factor in x at end-of-timestep (1=no clip, 0=full clip)', 'nondim')
      CS%id_phi_y_DG = register_diag_field('ice_shelf_model','phi_y_DG',CS%diag%axesT1, Time, &
         'DG(1) slope-limiter factor in y at end-of-timestep (1=no clip, 0=full clip)', 'nondim')
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
        call calc_shelf_driving_stress_DG(CS, ISS, G, US, CS%taudx_shelf, CS%taudy_shelf, CS%OD_av)
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
    if (CS%id_h_x > 0) call post_data(CS%id_h_x, CS%h_x, CS%diag)
    if (CS%id_h_y > 0) call post_data(CS%id_h_y, CS%h_y, CS%diag)
    if (CS%id_phi_x_DG > 0 .and. associated(CS%phi_x_DG)) &
        call post_data(CS%id_phi_x_DG, CS%phi_x_DG, CS%diag)
    if (CS%id_phi_y_DG > 0 .and. associated(CS%phi_y_DG)) &
        call post_data(CS%id_phi_y_DG, CS%phi_y_DG, CS%diag)
    if (CS%id_phi_x_FV > 0 .and. associated(CS%phi_x_FV)) &
        call post_data(CS%id_phi_x_FV, CS%phi_x_FV, CS%diag)
    if (CS%id_phi_y_FV > 0 .and. associated(CS%phi_y_FV)) &
        call post_data(CS%id_phi_y_FV, CS%phi_y_FV, CS%diag)
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
    ! DG(1) unsplit advection with SSP-RK2
    call ice_shelf_advect_DG1(CS, ISS, G, time_step, ISS%hmask, ISS%h_shelf, CS%h_x, CS%h_y, uh_ice, vh_ice)
    ! call pass_var(ISS%h_shelf, G%domain)
    ! call pass_var(CS%h_x, G%domain)
    ! call pass_var(CS%h_y, G%domain)
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
      call ice_shelf_min_thickness_calve(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, &
                                         CS%min_thickness_simple_calve)
    endif
    if (CS%calve_to_mask) then
      call calve_to_mask(G, ISS%h_shelf, ISS%area_shelf_h, ISS%hmask, CS%calve_mask)
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
  call compute_ground_frac(CS, ISS, G, H_node)

  ! Calculate RHS
  if (CS%use_DG_thickness) then
    call calc_shelf_driving_stress_DG(CS, ISS, G, US, taudx, taudy, CS%OD_av)
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
subroutine ice_shelf_min_thickness_calve(G, h_shelf, area_shelf_h, hmask, thickness_calve, halo)
  type(ocean_grid_type), intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real,                  intent(in)    :: thickness_calve !< The thickness at which to trigger calving [Z ~> m].
  integer,     optional, intent(in)    :: halo  !< The number of halo points to use.  If not present,
                                                !! work on the entire data domain.
  integer :: i, j, is, ie, js, je

  if (present(halo)) then
    is = G%isc - halo ; ie = G%iec + halo ; js = G%jsc - halo ; je = G%jec + halo
  else
    is = G%isd ; ie = G%ied ; js = G%jsd ; je = G%jed
  endif

  do j=js,je ; do i=is,ie
!    if ((h_shelf(i,j) < CS%thickness_calve) .and. (hmask(i,j) == 1) .and. &
!        (CS%ground_frac(i,j) == 0.0)) then
    if ((h_shelf(i,j) < thickness_calve) .and. (area_shelf_h(i,j) > 0.)) then
      h_shelf(i,j) = 0.0
      area_shelf_h(i,j) = 0.0
      hmask(i,j) = 0.0
    endif
  enddo ; enddo

end subroutine ice_shelf_min_thickness_calve

subroutine calve_to_mask(G, h_shelf, area_shelf_h, hmask, calve_mask)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: area_shelf_h !< The area per cell covered by
                                                             !! the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), intent(inout) :: hmask !< A mask indicating which tracer points are
                                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in)    :: calve_mask !< A mask that indicates where the ice
                                                             !! shelf can exist, and where it will calve.

  integer                        :: i,j

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((calve_mask(i,j) == 0.0) .and. (hmask(i,j) /= 0.0)) then
      h_shelf(i,j) = 0.0
      area_shelf_h(i,j) = 0.0
      hmask(i,j) = 0.0
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
          neumann_val = ((.5 * grav) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2 - rhow * CS%bed_elev(i,j)**2))
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
          ! DG mode: per-Gauss-point grounding check and fB computation
          if (do_DG) then
            h_gp = max(h_shelf(i,j) + ((CS%h_x(i,j)*(xquad(iq)-0.5)) + (CS%h_y(i,j)*(xquad(jq)-0.5))), &
                       CS%min_h_shelf)
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
          call CG_action_subgrid_basal(CS, G, US, Phisub, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
              bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
              use_DG=.true., h_shelf_cell=h_shelf(i,j), &
              h_x_cell=CS%h_x(i,j), h_y_cell=CS%h_y(i,j), &
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
                                   use_DG, h_shelf_cell, h_x_cell, h_y_cell, bed_corners)
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
  real,          optional, intent(in) :: h_x_cell      !< DG x-slope moment [Z ~> m]
  real,          optional, intent(in) :: h_y_cell      !< DG y-slope moment [Z ~> m]
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
        hloc = max(h_shelf_cell + ((h_x_cell*xi_sub) + (h_y_cell*eta_sub)), CS%min_h_shelf)
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
          h_gp = max(h_shelf(i,j) + ((CS%h_x(i,j)*(xquad(iq)-0.5)) + (CS%h_y(i,j)*(xquad(jq)-0.5))), &
                     CS%min_h_shelf)
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
        call CG_diagonal_subgrid_basal(CS, G, US, Phisub, Hcell, &
            u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
            CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
            G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
            use_DG=.true., h_shelf_cell=h_shelf(i,j), &
            h_x_cell=CS%h_x(i,j), h_y_cell=CS%h_y(i,j), &
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
                                     use_DG, h_shelf_cell, h_x_cell, h_y_cell, bed_corners)
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
  real,          optional, intent(in) :: h_x_cell      !< DG x-slope moment [Z ~> m]
  real,          optional, intent(in) :: h_y_cell      !< DG y-slope moment [Z ~> m]
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
        ! xi_sub = a_right(qx,i) - 0.5; marginal sum of Phisub over the l index gives a_right(qx,i).
        xi_sub  = (Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)) - 0.5
        eta_sub = (Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)) - 0.5
        hloc = max(h_shelf_cell + ((h_x_cell*xi_sub) + (h_y_cell*eta_sub)), CS%min_h_shelf)
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
    call interpolate_H_to_B(G, ISS%h_shelf, ISS%hmask, H_node, CS%min_h_shelf)

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
  real :: xi_gp, eta_gp ! Reference element coordinates at Gauss point [nondim]
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
            xi_gp  = xquad(iq) - 0.5  ! map [0,1] Gauss node to [-0.5,0.5] ref coords
            eta_gp = xquad(jq) - 0.5
            h_gp = max(ISS%h_shelf(i,j) + ((CS%h_x(i,j)*xi_gp) + (CS%h_y(i,j)*eta_gp)), &
                       CS%min_h_shelf)
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
  real :: xi_sub, eta_sub ! DG reference coords at sub-qp ([-0.5,0.5]) [nondim]
  real :: h_ip, bed_ip   ! Ice thickness and bed elevation at a sub-IP [Z ~> m]
  real :: bed_corners(2,2) ! Bed elevation at the 4 B-grid corners of a cell [Z ~> m]
  real :: H_corners(2,2)   ! Ice thickness at the 4 B-grid corners (non-DG path) [Z ~> m]
  integer :: i, j, isub, jsub, iq, jq, n_total, n_grounded
  integer :: isc, iec, jsc, jec

  if (.not. CS%GL_regularize) return

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
          xi_sub  = (CS%Phisub(iq,jq,isub,jsub,2,1) + CS%Phisub(iq,jq,isub,jsub,2,2)) - 0.5
          eta_sub = (CS%Phisub(iq,jq,isub,jsub,1,2) + CS%Phisub(iq,jq,isub,jsub,2,2)) - 0.5
          h_ip = max(ISS%h_shelf(i,j) + ((CS%h_x(i,j)*xi_sub) + (CS%h_y(i,j)*eta_sub)), &
                     CS%min_h_shelf)
        else
          h_ip = ((CS%Phisub(iq,jq,isub,jsub,1,1)*H_corners(1,1)) + (CS%Phisub(iq,jq,isub,jsub,2,2)*H_corners(2,2))) + &
                 ((CS%Phisub(iq,jq,isub,jsub,2,1)*H_corners(2,1)) + (CS%Phisub(iq,jq,isub,jsub,1,2)*H_corners(1,2)))
          h_ip = max(h_ip, CS%min_h_shelf)
        endif

        bed_ip = ((CS%Phisub(iq,jq,isub,jsub,1,1)*bed_corners(1,1)) + (CS%Phisub(iq,jq,isub,jsub,2,2)*bed_corners(2,2))) + &
                 ((CS%Phisub(iq,jq,isub,jsub,2,1)*bed_corners(2,1)) + (CS%Phisub(iq,jq,isub,jsub,1,2)*bed_corners(1,2)))

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
    if (J>1) then
      a = (G%dxCv(i,J-1) * yquad_m(qpoint)) + (G%dxCv(i,J) * yquad(qpoint)) ! d(x)/d(x*)
    else
      a = G%dxCv(i,J) !* yquad(qpoint) ! d(x)/d(x*)
    endif
    if (I>1) then
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

    ! d(x)/d(x*)
    if (J>1) then
      a = 0.5 * (G%dxCv(i,J-1) + G%dxCv(i,J))
    else
      a = G%dxCv(i,J)
    endif

    ! d(y)/d(y*)
    if (I>1) then
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
  deallocate(CS%h_x, CS%h_y)
  if (associated(CS%phi_x_DG)) deallocate(CS%phi_x_DG)
  if (associated(CS%phi_y_DG)) deallocate(CS%phi_y_DG)
  if (associated(CS%phi_x_FV)) deallocate(CS%phi_x_FV)
  if (associated(CS%phi_y_FV)) deallocate(CS%phi_y_FV)
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

!> Advect ice shelf thickness using a DG(1) (discontinuous Galerkin, linear polynomial)
!! method. Each cell stores 3 DOFs: cell-average h_shelf plus slope moments h_x, h_y.
!! The polynomial in cell (i,j) is: h(xi,eta) = h_shelf(i,j) + h_x(i,j)*xi + h_y(i,j)*eta
!! where (xi,eta) in [-0.5,0.5] are local cell coordinates.
!! Uses SSP-RK2 time stepping with upwind numerical fluxes and a minmod slope limiter.
!! This is an unsplit 2D scheme (no directional splitting).
subroutine ice_shelf_advect_DG1(CS, ISS, G, time_step, hmask, h_shelf, h_x, h_y, uh_ice, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real,                   intent(in)    :: time_step !< The time step for this update [T ~> s]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_shelf !< The cell-averaged ice shelf thickness [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_x !< DG(1) x-slope moment [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_y !< DG(1) y-slope moment [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< The accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< The accumulated meridional ice volume flux [Z L2 ~> m3]

  ! Local variables for SSP-RK2 stages
  real, dimension(SZDI_(G),SZDJ_(G)) :: h0, hx0, hy0       ! Initial state [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: h1, hx1, hy1       ! After RK stage 1 [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: Rhs_h, Rhs_hx, Rhs_hy ! Spatial operator [Z T-1 ~> m s-1]
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  uh_ice(:,:) = 0.0
  vh_ice(:,:) = 0.0

  ! Apply boundary values
  do j=jsd,jed ; do i=isd,ied
    if (CS%h_bdry_val(i,j) /= 0.0) then
      h_shelf(i,j) = CS%h_bdry_val(i,j)
      h_x(i,j) = 0.0 ; h_y(i,j) = 0.0
    endif
  enddo ; enddo

  ! Save initial state
  do j=jsd,jed ; do i=isd,ied
    h0(i,j) = h_shelf(i,j) ; hx0(i,j) = h_x(i,j) ; hy0(i,j) = h_y(i,j)
  enddo ; enddo
  call pass_var(h0, G%domain)

  ! --- SSP-RK2 Stage 1: h1 = h0 + dt * L(h0) ---
  call DG1_slope_limit(G, h0, hx0, hy0, hmask, CS%h_bdry_val, &
                       CS%dg1_limiter_choice, CS%dg1_limiter_M)
  call pass_vector(hx0, hy0, G%domain, TO_ALL, AGRID)
  call DG1_spatial_operator(CS, G, hmask, h0, hx0, hy0, Rhs_h, Rhs_hx, Rhs_hy, uh_ice, vh_ice)

  do j=jsc,jec ; do i=isc,iec
    if (hmask(i,j) == 1) then
      h1(i,j)  = h0(i,j)  + time_step * Rhs_h(i,j)
      hx1(i,j) = hx0(i,j) + time_step * Rhs_hx(i,j)
      hy1(i,j) = hy0(i,j) + time_step * Rhs_hy(i,j)
    else
      h1(i,j) = h0(i,j) ; hx1(i,j) = hx0(i,j) ; hy1(i,j) = hy0(i,j)
    endif
  enddo ; enddo
  call pass_var(h1, G%domain)
  call pass_vector(hx1, hy1, G%domain, TO_ALL, AGRID)

  ! --- SSP-RK2 Stage 2: h_new = 0.5*h0 + 0.5*(h1 + dt * L(h1)) ---
  call DG1_slope_limit(G, h1, hx1, hy1, hmask, CS%h_bdry_val, &
                       CS%dg1_limiter_choice, CS%dg1_limiter_M)
  call DG1_spatial_operator(CS, G, hmask, h1, hx1, hy1, Rhs_h, Rhs_hx, Rhs_hy, uh_ice, vh_ice)

  do j=jsc,jec ; do i=isc,iec
    if (hmask(i,j) == 1) then
      h_shelf(i,j) = 0.5 * h0(i,j) + 0.5 * (h1(i,j)  + time_step * Rhs_h(i,j))
      h_x(i,j)     = 0.5 * hx0(i,j)+ 0.5 * (hx1(i,j) + time_step * Rhs_hx(i,j))
      h_y(i,j)     = 0.5 * hy0(i,j)+ 0.5 * (hy1(i,j) + time_step * Rhs_hy(i,j))
    endif
  enddo ; enddo
  call pass_var(h_shelf, G%domain)
  call pass_vector(h_x, h_y, G%domain, TO_ALL, AGRID)

  ! Final slope limit. Capture the limiter factors here so they reflect the
  ! limiter's effect on the state stored at end-of-timestep.
  call DG1_slope_limit(G, h_shelf, h_x, h_y, hmask, CS%h_bdry_val, &
                       CS%dg1_limiter_choice, CS%dg1_limiter_M, &
                       phi_x_out=CS%phi_x_DG, phi_y_out=CS%phi_y_DG)

  call pass_var(h_shelf, G%domain)
  call pass_vector(h_x, h_y, G%domain, TO_ALL, AGRID)

  ! Scale uh_ice, vh_ice: the spatial operator accumulated fluxes from both RK stages,
  ! so average them (SSP-RK2 gives equal weight to each stage)
  do j=jsc,jec ; do I=isc-1,iec
    uh_ice(I,j) = 0.5 * uh_ice(I,j)
  enddo ; enddo
  do J=jsc-1,jec ; do i=isc,iec
    vh_ice(i,J) = 0.5 * vh_ice(i,J)
  enddo ; enddo

  call pass_vector(uh_ice, vh_ice, G%domain, TO_ALL, CGRID_NE)

end subroutine ice_shelf_advect_DG1


!> Compute the DG(1) spatial operator (right-hand side) for the thickness evolution
!! equation dh/dt = -div(u*h). Computes the flux divergence contribution to each
!! moment (cell average, x-slope, y-slope) using upwind numerical fluxes with
!! 2-point Gauss quadrature along each face.
subroutine DG1_spatial_operator(CS, G, hmask, h_bar, h_x, h_y, Rhs_h, Rhs_hx, Rhs_hy, uh_ice, vh_ice)
  type(ice_shelf_dyn_CS), intent(in)    :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf mask
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_bar !< Cell-averaged thickness [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_x  !< DG x-slope moment [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_y  !< DG y-slope moment [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(out)   :: Rhs_h  !< RHS for cell average [Z T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(out)   :: Rhs_hx !< RHS for x-slope moment [Z T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(out)   :: Rhs_hy !< RHS for y-slope moment [Z T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJ_(G)), &
                          intent(inout) :: uh_ice !< Accumulated zonal ice volume flux [Z L2 ~> m3]
  real, dimension(SZDI_(G),SZDJB_(G)), &
                          intent(inout) :: vh_ice !< Accumulated meridional ice volume flux [Z L2 ~> m3]

  ! Gauss quadrature points on [-0.5, 0.5] (2-point rule)
  real, parameter :: gp1 = -0.5/sqrt(3.0)  ! ~-0.2887 [nondim]
  real, parameter :: gp2 =  0.5/sqrt(3.0)  ! ~+0.2887 [nondim]
  real, parameter :: gw  =  0.5            ! Weight for each point (sums to 1) [nondim]

  ! Local variables
  real :: u_face        ! Normal velocity at a u-face [L T-1 ~> m s-1]
  real :: v_face        ! Normal velocity at a v-face [L T-1 ~> m s-1]
  real :: h_upwind      ! Upwind thickness at a Gauss point [Z ~> m]
  real :: flux_h        ! Flux contribution to cell average [Z L2 T-1 ~> m3 s-1]
  real :: flux_hx       ! Flux contribution to x-moment [Z L2 T-1 ~> m3 s-1]
  real :: flux_hy       ! Flux contribution to y-moment [Z L2 T-1 ~> m3 s-1]
  real :: eta_gp        ! Gauss point coordinate along face [nondim]
  real :: face_flux_total ! Total volume flux through a face [Z L2 T-1 ~> m3 s-1]
  real :: u_c, u_xi, u_eta ! Bilinear u modes on a cell from B-grid corners [L T-1 ~> m s-1]
  real :: v_c, v_xi, v_eta ! Bilinear v modes on a cell from B-grid corners [L T-1 ~> m s-1]
  ! Per-cell accumulators split by face orientation, so that the final reduction
  ! Rhs = Rhs_u + Rhs_v is a single binary add (order-independent under 90 deg
  ! rotation, which swaps the u-face and v-face contributions).
  real, dimension(SZDI_(G),SZDJ_(G)) :: Rhs_h_u, Rhs_h_v   ! Cell-mean RHS from u/v faces
  real, dimension(SZDI_(G),SZDJ_(G)) :: Rhs_hx_u, Rhs_hx_v ! x-moment RHS from u/v faces
  real, dimension(SZDI_(G),SZDJ_(G)) :: Rhs_hy_u, Rhs_hy_v ! y-moment RHS from u/v faces
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed, gp

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  Rhs_h(:,:) = 0.0 ; Rhs_hx(:,:) = 0.0 ; Rhs_hy(:,:) = 0.0
  Rhs_h_u(:,:)  = 0.0 ; Rhs_h_v(:,:)  = 0.0
  Rhs_hx_u(:,:) = 0.0 ; Rhs_hx_v(:,:) = 0.0
  Rhs_hy_u(:,:) = 0.0 ; Rhs_hy_v(:,:) = 0.0

  ! --- Zonal (east) face fluxes at I-faces ---
  ! Face I between cells (i,j) [left] and (i+1,j) [right]
  ! On this face, xi = +0.5 for the left cell, xi = -0.5 for the right cell
  ! eta varies along the face; Gauss points at eta = gp1, gp2
  do j=jsc,jec ; do I=isc-1,iec
    if (CS%u_face_mask(I,j) == 4.) then
      ! Specified flux boundary condition
      face_flux_total = G%dyCu(I,j) * CS%u_flux_bdry_val(I,j)
      uh_ice(I,j) = uh_ice(I,j) + face_flux_total
      ! For specified flux, add to cell average RHS only (no moment info in BC)
      if (hmask(i,j) == 1) &
        Rhs_h_u(i,j) = Rhs_h_u(i,j) - face_flux_total * G%IareaT(i,j)
      if (hmask(i+1,j) == 1) &
        Rhs_h_u(i+1,j) = Rhs_h_u(i+1,j) + face_flux_total * G%IareaT(i+1,j)
    elseif ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .or. &
            (hmask(i+1,j) == 1 .or. hmask(i+1,j) == 3)) then
      ! Normal velocity at this u-face (average of B-grid velocities)
      u_face = 0.5 * (CS%u_shelf(I,J-1) + CS%u_shelf(I,J))

      flux_h = 0.0 ; flux_hx = 0.0 ; flux_hy = 0.0

      ! 2-point Gauss quadrature along the face (in eta direction)
      do gp=1,2
        if (gp == 1) then ; eta_gp = gp1 ; else ; eta_gp = gp2 ; endif

        ! Upwind thickness at this Gauss point
        if (u_face > 0.0) then
          if (hmask(i,j) == 3) then
            h_upwind = CS%h_bdry_val(i,j)
          elseif (hmask(i,j) == 1) then
            ! Evaluate left cell polynomial at (xi=+0.5, eta=eta_gp)
            h_upwind = h_bar(i,j) + ((h_x(i,j) * 0.5) + (h_y(i,j) * eta_gp))
          else
            h_upwind = 0.0
          endif
        else
          if (hmask(i+1,j) == 3) then
            h_upwind = CS%h_bdry_val(i+1,j)
          elseif (hmask(i+1,j) == 1) then
            ! Evaluate right cell polynomial at (xi=-0.5, eta=eta_gp)
            h_upwind = h_bar(i+1,j) + ((h_x(i+1,j) * (-0.5)) + (h_y(i+1,j) * eta_gp))
          else
            h_upwind = 0.0
          endif
        endif

        h_upwind = max(h_upwind, 0.0)

        ! Numerical flux at this Gauss point: u * h_upwind * face_length
        ! Weight gw accounts for the Gauss quadrature (each point gets weight 0.5)
        flux_h  = flux_h  + gw * u_face * h_upwind * G%dyCu(I,j)
        ! Moment flux: weighted by the test function value at the face
        ! For the x-moment: test function xi = +0.5 at right face of left cell, -0.5 at left face of right cell
        ! For the y-moment: test function eta = eta_gp
        flux_hx = flux_hx + gw * u_face * h_upwind * G%dyCu(I,j) * 0.5   ! xi at face = +/-0.5
        flux_hy = flux_hy + gw * u_face * h_upwind * G%dyCu(I,j) * eta_gp
      enddo

      uh_ice(I,j) = uh_ice(I,j) + flux_h

      ! Subtract outgoing flux from left cell (i,j), add incoming flux to right cell (i+1,j)
      ! For cell average: Rhs_h -= flux / area  (flux leaving through east face)
      if (hmask(i,j) == 1) then
        Rhs_h_u(i,j)  = Rhs_h_u(i,j)  - flux_h * G%IareaT(i,j)
        ! For x-moment: the face value of the x test function is +0.5
        Rhs_hx_u(i,j) = Rhs_hx_u(i,j) - flux_hx * G%IareaT(i,j)
        ! For y-moment: weighted by eta_gp (already in flux_hy)
        Rhs_hy_u(i,j) = Rhs_hy_u(i,j) - flux_hy * G%IareaT(i,j)
      endif
      if (hmask(i+1,j) == 1) then
        Rhs_h_u(i+1,j)  = Rhs_h_u(i+1,j)  + flux_h * G%IareaT(i+1,j)
        ! For the right cell, xi at its left face = -0.5
        Rhs_hx_u(i+1,j) = Rhs_hx_u(i+1,j) + flux_hx * (-1.0) * G%IareaT(i+1,j)
        Rhs_hy_u(i+1,j) = Rhs_hy_u(i+1,j) + flux_hy * G%IareaT(i+1,j)
      endif
    endif
  enddo ; enddo

  ! --- Meridional (north) face fluxes at J-faces ---
  ! Face J between cells (i,j) [south] and (i,j+1) [north]
  ! On this face, eta = +0.5 for the south cell, eta = -0.5 for the north cell
  ! xi varies along the face; Gauss points at xi = gp1, gp2
  do J=jsc-1,jec ; do i=isc,iec
    if (CS%v_face_mask(i,J) == 4.) then
      ! Specified flux boundary condition
      face_flux_total = G%dxCv(i,J) * CS%v_flux_bdry_val(i,J)
      vh_ice(i,J) = vh_ice(i,J) + face_flux_total
      if (hmask(i,j) == 1) &
        Rhs_h_v(i,j) = Rhs_h_v(i,j) - face_flux_total * G%IareaT(i,j)
      if (hmask(i,j+1) == 1) &
        Rhs_h_v(i,j+1) = Rhs_h_v(i,j+1) + face_flux_total * G%IareaT(i,j+1)
    elseif ((hmask(i,j) == 1 .or. hmask(i,j) == 3) .or. &
            (hmask(i,j+1) == 1 .or. hmask(i,j+1) == 3)) then
      ! Normal velocity at this v-face
      v_face = 0.5 * (CS%v_shelf(I-1,J) + CS%v_shelf(I,J))

      flux_h = 0.0 ; flux_hx = 0.0 ; flux_hy = 0.0

      do gp=1,2
        if (gp == 1) then ; eta_gp = gp1 ; else ; eta_gp = gp2 ; endif
        ! Here eta_gp is used as the xi coordinate along the face

        if (v_face > 0.0) then
          if (hmask(i,j) == 3) then
            h_upwind = CS%h_bdry_val(i,j)
          elseif (hmask(i,j) == 1) then
            ! Evaluate south cell polynomial at (xi=eta_gp, eta=+0.5)
            h_upwind = h_bar(i,j) + ((h_x(i,j) * eta_gp) + (h_y(i,j) * 0.5))
          else
            h_upwind = 0.0
          endif
        else
          if (hmask(i,j+1) == 3) then
            h_upwind = CS%h_bdry_val(i,j+1)
          elseif (hmask(i,j+1) == 1) then
            ! Evaluate north cell polynomial at (xi=eta_gp, eta=-0.5)
            h_upwind = h_bar(i,j+1) + ((h_x(i,j+1) * eta_gp) + (h_y(i,j+1) * (-0.5)))
          else
            h_upwind = 0.0
          endif
        endif

        h_upwind = max(h_upwind, 0.0)

        flux_h  = flux_h  + gw * v_face * h_upwind * G%dxCv(i,J)
        flux_hx = flux_hx + gw * v_face * h_upwind * G%dxCv(i,J) * eta_gp  ! xi along face
        flux_hy = flux_hy + gw * v_face * h_upwind * G%dxCv(i,J) * 0.5     ! eta at face = +/-0.5
      enddo

      vh_ice(i,J) = vh_ice(i,J) + flux_h

      if (hmask(i,j) == 1) then
        Rhs_h_v(i,j)  = Rhs_h_v(i,j)  - flux_h * G%IareaT(i,j)
        Rhs_hx_v(i,j) = Rhs_hx_v(i,j) - flux_hx * G%IareaT(i,j)
        ! For south cell, eta at its north face = +0.5
        Rhs_hy_v(i,j) = Rhs_hy_v(i,j) - flux_hy * G%IareaT(i,j)
      endif
      if (hmask(i,j+1) == 1) then
        Rhs_h_v(i,j+1)  = Rhs_h_v(i,j+1)  + flux_h * G%IareaT(i,j+1)
        Rhs_hx_v(i,j+1) = Rhs_hx_v(i,j+1) + flux_hx * G%IareaT(i,j+1)
        ! For north cell, eta at its south face = -0.5
        Rhs_hy_v(i,j+1) = Rhs_hy_v(i,j+1) + flux_hy * (-1.0) * G%IareaT(i,j+1)
      endif
    endif
  enddo ; enddo

  ! Volume integral term for the slope moments:
  ! The DG(1) weak form includes a volume integral: integral(u*h * d(test)/dx, dA)
  ! For the x-moment with test function xi: d(xi)/dx = 1/dx_cell, so the volume term is
  !   (1/A) * integral(u*h * (1/dx_cell), dA) which approximated at cell center gives:
  !   u_center * h_bar / dx_cell  (but this is already captured by the face fluxes in the
  !   DG formulation after integration by parts). The face flux terms above already include
  !   the complete weak-form contributions including the volume integral, since we used
  !   integration by parts: integral(div(u*h)*test, dA) = boundary(u*h*test, ds) - integral(u*h*grad(test), dA)
  ! The above face flux terms are the boundary integral. We need to add the volume integral.

  ! Reduce the per-face accumulators into the destination Rhs arrays. The single
  ! binary add Rhs = Rhs_u + Rhs_v is FP-commutative, so the result is bit-identical
  ! under a 90 deg rotation (which swaps the u-face and v-face contributions).
  do j=jsc,jec ; do i=isc,iec
    if (hmask(i,j) == 1) then
      Rhs_h(i,j)  = Rhs_h_u(i,j)  + Rhs_h_v(i,j)
      Rhs_hx(i,j) = Rhs_hx_u(i,j) + Rhs_hx_v(i,j)
      Rhs_hy(i,j) = Rhs_hy_u(i,j) + Rhs_hy_v(i,j)
    endif
  enddo ; enddo

  ! Volume integral contribution: + integral(u*h * d(phi)/dx, dA) for each test function.
  ! For test xi: d(xi)/dx = 1/dx, d(xi)/dy = 0; for test eta: vice versa.
  ! With bilinear u(xi,eta) = u_c + u_xi*xi + u_eta*eta + u_xieta*xi*eta from B-grid
  ! corners, and linear h(xi,eta) = h_bar + h_x*xi + h_y*eta, the integral over
  ! [-0.5,0.5]^2 evaluates to:
  !   integral(u*h, dxi deta) = u_c*h_bar + (u_xi*h_x + u_eta*h_y)/12
  ! (odd-power terms vanish; the xi*eta cross term contributes nothing because it
  ! pairs only with odd-power h modes). Same form for v in the y-moment.
  do j=jsc,jec ; do i=isc,iec
    if (hmask(i,j) == 1) then
      ! Diagonal + off-diagonal grouping so the 4-corner mean is invariant under
      ! 90 deg rotation (which permutes corners within these two pairs).
      u_c   = 0.25 * ((CS%u_shelf(I-1,J-1) + CS%u_shelf(I,J)) + &
                      (CS%u_shelf(I,J-1)   + CS%u_shelf(I-1,J)))
      u_xi  = 0.5  * ((CS%u_shelf(I,J-1)   + CS%u_shelf(I,J)) - &
                      (CS%u_shelf(I-1,J-1) + CS%u_shelf(I-1,J)))
      u_eta = 0.5  * ((CS%u_shelf(I-1,J)   + CS%u_shelf(I,J)) - &
                      (CS%u_shelf(I-1,J-1) + CS%u_shelf(I,J-1)))
      v_c   = 0.25 * ((CS%v_shelf(I-1,J-1) + CS%v_shelf(I,J)) + &
                      (CS%v_shelf(I,J-1)   + CS%v_shelf(I-1,J)))
      v_xi  = 0.5  * ((CS%v_shelf(I,J-1)   + CS%v_shelf(I,J)) - &
                      (CS%v_shelf(I-1,J-1) + CS%v_shelf(I-1,J)))
      v_eta = 0.5  * ((CS%v_shelf(I-1,J)   + CS%v_shelf(I,J)) - &
                      (CS%v_shelf(I-1,J-1) + CS%v_shelf(I,J-1)))

      Rhs_hx(i,j) = Rhs_hx(i,j) + &
        ((u_c * h_bar(i,j)) + (((u_xi * h_x(i,j)) + (u_eta * h_y(i,j))) / 12.0)) * G%IdxT(i,j)
      Rhs_hy(i,j) = Rhs_hy(i,j) + &
        ((v_c * h_bar(i,j)) + (((v_xi * h_x(i,j)) + (v_eta * h_y(i,j))) / 12.0)) * G%IdyT(i,j)
    endif
  enddo ; enddo

  ! Scale moment RHS by 12 for the inverse mass matrix of linear DG basis on [-0.5,0.5]:
  ! The mass matrix for basis {1, xi, eta} on [-0.5,0.5]^2 is diag(1, 1/12, 1/12)
  ! So the inverse mass matrix scales the slope moment RHS by 12.
  do j=jsc,jec ; do i=isc,iec
    if (hmask(i,j) == 1) then
      Rhs_hx(i,j) = 12.0 * Rhs_hx(i,j)
      Rhs_hy(i,j) = 12.0 * Rhs_hy(i,j)
    endif
  enddo ; enddo

end subroutine DG1_spatial_operator


!> Apply a slope limiter to the DG(1) slope moments to control oscillations.
!! Operates per-direction on h_x and h_y so the algorithm is symmetric in
!! x and y (rotation-invariant under 90-degree grid rotations).
!! Dispatches on limiter_choice:
!!  0 = no limiting,
!!  1 = minmod (legacy),
!!  2 = Venkatakrishnan with TVB-Cockburn-Shu eps band eps = M*dx**2.
!! When phi_x_out / phi_y_out are present, the limiter factors that were
!! applied to h_x, h_y are returned (1 = no clipping, 0 = full clip).
subroutine DG1_slope_limit(G, h_bar, h_x, h_y, hmask, h_bdry_val, &
                           limiter_choice, vk_M, phi_x_out, phi_y_out)
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_bar !< Cell-averaged thickness [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_x  !< DG x-slope moment [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_y  !< DG y-slope moment [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: hmask !< Ice shelf mask
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(in)    :: h_bdry_val !< Dirichlet boundary thickness, used as
                                              !! a face value at hmask==3 cells [Z ~> m]
  integer,                intent(in)    :: limiter_choice !< 0=none, 1=minmod, 2=Venkatakrishnan
  real,                   intent(in)    :: vk_M !< TVB curvature bound [Z L-2 ~> m-1]
                                              !! used when limiter_choice == 2.
  real, dimension(SZDI_(G),SZDJ_(G)), optional, &
                          intent(out)   :: phi_x_out !< x limiter factor [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), optional, &
                          intent(out)   :: phi_y_out !< y limiter factor [nondim]

  real :: diff_E, diff_W, diff_N, diff_S ! Neighbor differences [Z ~> m]
  real :: h_E_eff, h_W_eff, h_N_eff, h_S_eff ! Effective neighbour means [Z ~> m]
  real :: h_max_x, h_min_x, h_max_y, h_min_y ! Cell-neighbourhood envelope [Z ~> m]
  real :: delta1_E, delta1_W, delta1_N, delta1_S ! Allowed face increments [Z ~> m]
  real :: delta2_E, delta2_W, delta2_N, delta2_S ! Predicted face increments [Z ~> m]
  real :: phi_E, phi_W, phi_N, phi_S, phi_x, phi_y ! Per-face/cell Venk factors [nondim]
  real :: eps_sq_x, eps_sq_y             ! Smooth-extremum band squared [Z2 ~> m2]
  integer :: i, j, isc, iec, jsc, jec
  logical :: valid_E, valid_W, valid_N, valid_S

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  if (present(phi_x_out)) phi_x_out(:,:) = 1.0
  if (present(phi_y_out)) phi_y_out(:,:) = 1.0

  if (limiter_choice == 0) return ! "none"

  do j=jsc,jec ; do i=isc,iec
    if (hmask(i,j) /= 1) then
      h_x(i,j) = 0.0 ; h_y(i,j) = 0.0
      ! if (present(phi_x_out)) phi_x_out(i,j) = 0.0
      ! if (present(phi_y_out)) phi_y_out(i,j) = 0.0
      cycle
    endif

    valid_E = (hmask(i+1,j) == 1 .or. hmask(i+1,j) == 3)
    valid_W = (hmask(i-1,j) == 1 .or. hmask(i-1,j) == 3)
    valid_N = (hmask(i,j+1) == 1 .or. hmask(i,j+1) == 3)
    valid_S = (hmask(i,j-1) == 1 .or. hmask(i,j-1) == 3)

    if (limiter_choice == 1) then
      ! ---- Minmod (legacy) ----
      ! At hmask==3 neighbours, h_bdry_val is the FACE value at distance dx/2
      ! (see ice_shelf_advect_thickness_x). The cell-mean difference across one
      ! cell width is therefore 2*(h_bar(i) - h_bdry_val) on that side.

      ! X-direction limiter. Use a one-sided difference when only one neighbour
      ! is valid (e.g. cell adjacent to the ice front) so the slope is preserved
      ! rather than zeroed: zeroing biases the cell-face value used for outflow
      ! and pollutes the SSA driving stress at the adjacent corner nodes.
      if (valid_E) then
        if (hmask(i+1,j) == 3) then
          diff_E = 2.0 * (h_bdry_val(i+1,j) - h_bar(i,j))
        else
          diff_E = h_bar(i+1,j) - h_bar(i,j)
        endif
      endif
      if (valid_W) then
        if (hmask(i-1,j) == 3) then
          diff_W = 2.0 * (h_bar(i,j) - h_bdry_val(i-1,j))
        else
          diff_W = h_bar(i,j) - h_bar(i-1,j)
        endif
      endif
      if (valid_E .and. valid_W) then
        h_x(i,j) = minmod3(h_x(i,j), diff_E, diff_W)
      elseif (valid_W) then
        h_x(i,j) = minmod2(h_x(i,j), diff_W)
      elseif (valid_E) then
        h_x(i,j) = minmod2(h_x(i,j), diff_E)
      else
        h_x(i,j) = 0.0
      endif

      ! Y-direction limiter (same one-sided handling as X)
      if (valid_N) then
        if (hmask(i,j+1) == 3) then
          diff_N = 2.0 * (h_bdry_val(i,j+1) - h_bar(i,j))
        else
          diff_N = h_bar(i,j+1) - h_bar(i,j)
        endif
      endif
      if (valid_S) then
        if (hmask(i,j-1) == 3) then
          diff_S = 2.0 * (h_bar(i,j) - h_bdry_val(i,j-1))
        else
          diff_S = h_bar(i,j) - h_bar(i,j-1)
        endif
      endif
      if (valid_N .and. valid_S) then
        h_y(i,j) = minmod3(h_y(i,j), diff_N, diff_S)
      elseif (valid_S) then
        h_y(i,j) = minmod2(h_y(i,j), diff_S)
      elseif (valid_N) then
        h_y(i,j) = minmod2(h_y(i,j), diff_N)
      else
        h_y(i,j) = 0.0
      endif

    else
      ! ---- Venkatakrishnan with TVB-Cockburn-Shu eps band ----
      ! eps = M * dx_local**2 (separately in x and y). Per-direction structure
      ! is symmetric in x and y so the algorithm is invariant under 90-degree
      ! grid rotations.

      ! Effective neighbour cell-mean values via face-mirror at hmask==3.
      if (valid_E) then
        if (hmask(i+1,j) == 3) then
          h_E_eff = 2.0 * h_bdry_val(i+1,j) - h_bar(i,j)
        else
          h_E_eff = h_bar(i+1,j)
        endif
      endif
      if (valid_W) then
        if (hmask(i-1,j) == 3) then
          h_W_eff = 2.0 * h_bdry_val(i-1,j) - h_bar(i,j)
        else
          h_W_eff = h_bar(i-1,j)
        endif
      endif
      if (valid_N) then
        if (hmask(i,j+1) == 3) then
          h_N_eff = 2.0 * h_bdry_val(i,j+1) - h_bar(i,j)
        else
          h_N_eff = h_bar(i,j+1)
        endif
      endif
      if (valid_S) then
        if (hmask(i,j-1) == 3) then
          h_S_eff = 2.0 * h_bdry_val(i,j-1) - h_bar(i,j)
        else
          h_S_eff = h_bar(i,j-1)
        endif
      endif

      ! ---- X-direction limiter ----
      ! Cell-neighbourhood envelope (cell mean + valid x-neighbours):
      h_max_x = h_bar(i,j) ; h_min_x = h_bar(i,j)
      if (valid_E) then
        h_max_x = max(h_max_x, h_E_eff) ; h_min_x = min(h_min_x, h_E_eff)
      endif
      if (valid_W) then
        h_max_x = max(h_max_x, h_W_eff) ; h_min_x = min(h_min_x, h_W_eff)
      endif

      ! Predicted face increments (basis: face = h_bar +/- h_x/2):
      delta2_E = +0.5 * h_x(i,j)
      delta2_W = -0.5 * h_x(i,j)

      ! TVB-Cockburn-Shu eps band, cell-local in dx:
      eps_sq_x = (vk_M * G%dxT(i,j)*G%dxT(i,j))**2

      ! Per-face Venkatakrishnan factor:
      if (valid_E) then
        if (delta2_E > 0.0) then
          delta1_E = h_max_x - h_bar(i,j)
        else
          delta1_E = h_min_x - h_bar(i,j)
        endif
        phi_E = venk_factor(delta1_E, delta2_E, eps_sq_x)
      else
        phi_E = huge(1.0)
      endif
      if (valid_W) then
        if (delta2_W > 0.0) then
          delta1_W = h_max_x - h_bar(i,j)
        else
          delta1_W = h_min_x - h_bar(i,j)
        endif
        phi_W = venk_factor(delta1_W, delta2_W, eps_sq_x)
      else
        phi_W = huge(1.0)
      endif
      if (valid_E .or. valid_W) then
        phi_x = min(phi_E, phi_W)
      else
        phi_x = 0.0
      endif
      h_x(i,j) = phi_x * h_x(i,j)
      if (present(phi_x_out)) phi_x_out(i,j) = phi_x

      ! ---- Y-direction limiter ----
      ! Cell-neighbourhood envelope (cell mean + valid y-neighbours):
      h_max_y = h_bar(i,j) ; h_min_y = h_bar(i,j)
      if (valid_N) then
        h_max_y = max(h_max_y, h_N_eff) ; h_min_y = min(h_min_y, h_N_eff)
      endif
      if (valid_S) then
        h_max_y = max(h_max_y, h_S_eff) ; h_min_y = min(h_min_y, h_S_eff)
      endif

      delta2_N = +0.5 * h_y(i,j)
      delta2_S = -0.5 * h_y(i,j)

      eps_sq_y = (vk_M * G%dyT(i,j)*G%dyT(i,j))**2

      if (valid_N) then
        if (delta2_N > 0.0) then
          delta1_N = h_max_y - h_bar(i,j)
        else
          delta1_N = h_min_y - h_bar(i,j)
        endif
        phi_N = venk_factor(delta1_N, delta2_N, eps_sq_y)
      else
        phi_N = huge(1.0)
      endif
      if (valid_S) then
        if (delta2_S > 0.0) then
          delta1_S = h_max_y - h_bar(i,j)
        else
          delta1_S = h_min_y - h_bar(i,j)
        endif
        phi_S = venk_factor(delta1_S, delta2_S, eps_sq_y)
      else
        phi_S = huge(1.0)
      endif
      if (valid_N .or. valid_S) then
        phi_y = min(phi_N, phi_S)
      else
        phi_y = 0.0
      endif
      h_y(i,j) = phi_y * h_y(i,j)
      if (present(phi_y_out)) phi_y_out(i,j) = phi_y
    endif
  enddo ; enddo

end subroutine DG1_slope_limit


!> Smooth Venkatakrishnan limiter factor for one face of a DG cell.
!! Returns phi in [0, 1]: 1 means no clipping (predicted face value lies
!! safely inside the cell-neighbourhood envelope, or both delta1 and delta2
!! are small relative to eps so the face is treated as a smooth extremum);
!! 0 means full clip; intermediate values smoothly clip.
pure real function venk_factor(delta1, delta2, eps_sq)
  real, intent(in) :: delta1   !< Allowed increment to local max/min [Z ~> m]
  real, intent(in) :: delta2   !< Predicted increment from DG slope [Z ~> m]
  real, intent(in) :: eps_sq   !< Smooth-extremum band squared [Z2 ~> m2]
  real :: num, denom

  if (delta2 == 0.0) then
    venk_factor = 1.0
    return
  endif

  num   = (delta1*delta1 + eps_sq) * delta2 + 2.0 * delta2*delta2 * delta1
  denom = delta2 * (delta1*delta1 + 2.0*delta2*delta2 + delta1*delta2 + eps_sq)

  if (denom == 0.0) then
    venk_factor = 1.0
  else
    venk_factor = num / denom
  endif
  ! Clamp to [0, 1]; the analytic value is in this range when delta1 and
  ! delta2 agree in sign, but eps>0 can produce tiny floating-point drift.
  venk_factor = max(0.0, min(1.0, venk_factor))
end function venk_factor


!> Read DG(1) slope-limiter runtime parameters and map the choice string
!! to the integer enum stored in CS%dg1_limiter_choice.
subroutine read_dg1_limiter_params(param_file, mdl, CS, US)
  type(param_file_type),   intent(in)    :: param_file !< Parameter file
  character(len=*),        intent(in)    :: mdl        !< Module name
  type(ice_shelf_dyn_CS),  intent(inout) :: CS         !< Ice-shelf control structure
  type(unit_scale_type),   intent(in)    :: US         !< Unit scaling structure

  character(len=40) :: limiter_str

  call get_param(param_file, mdl, "DG1_LIMITER", limiter_str, &
                 "Slope limiter for DG(1) ice thickness. One of: "//&
                 "'venkatakrishnan' (default), 'minmod' (legacy, biased "//&
                 "on smooth concave-monotone flow), 'none' (no limiting). "//&
                 "Venkatakrishnan returns the DG-evolved slope unchanged when "//&
                 "the predicted face value lies safely inside the cell-"//&
                 "neighbourhood envelope, smoothly clips when it doesn't, "//&
                 "and protects smooth interior extrema via the eps band set "//&
                 "by DG1_LIMITER_M.", &
                 default="venkatakrishnan", do_not_log=.not.CS%use_DG_thickness)
  select case (trim(limiter_str))
  case ("none");            CS%dg1_limiter_choice = 0
  case ("minmod");          CS%dg1_limiter_choice = 1
  case ("venkatakrishnan"); CS%dg1_limiter_choice = 2
  case default
    call MOM_error(FATAL, "read_dg1_limiter_params: DG1_LIMITER must be "//&
                          "one of: none, minmod, venkatakrishnan.")
  end select

  call get_param(param_file, mdl, "DG1_LIMITER_M", CS%dg1_limiter_M, &
                 "TVB-style curvature bound for the Venkatakrishnan smooth-"//&
                 "extremum protection band. The local band is set as eps = "//&
                 "M * dx_local^2; M is the user's upper bound on |d^2 h / "//&
                 "dx^2| in smooth regions of the solution. Larger M widens "//&
                 "the bypass (more smooth extrema protected); smaller M "//&
                 "tightens it (more clipping). M=0 recovers pure smooth "//&
                 "Venkatakrishnan with no extremum protection; M very "//&
                 "large effectively disables the limiter. Only used when "//&
                 "DG1_LIMITER=venkatakrishnan.", &
                 units="m-1", default=1.0e-6, &
                 scale=US%m_to_Z/(US%m_to_L*US%m_to_L), &
                 do_not_log=(.not.CS%use_DG_thickness) .or. &
                            (CS%dg1_limiter_choice /= 2))

end subroutine read_dg1_limiter_params


!> Three-argument minmod function used by the DG slope limiter.
!! Returns 0 if the arguments have different signs, otherwise returns
!! the argument with the smallest absolute value.
pure real function minmod3(a, b, c)
  real, intent(in) :: a !< First argument [Z ~> m]
  real, intent(in) :: b !< Second argument [Z ~> m]
  real, intent(in) :: c !< Third argument [Z ~> m]

  if (a > 0.0 .and. b > 0.0 .and. c > 0.0) then
    minmod3 = min(a, b, c)
  elseif (a < 0.0 .and. b < 0.0 .and. c < 0.0) then
    minmod3 = max(a, b, c)  ! max of negatives = smallest magnitude
  else
    minmod3 = 0.0
  endif
end function minmod3


!> Two-argument minmod: returns 0 if arguments have different signs, otherwise
!! the argument with the smaller absolute value. Used for one-sided slope
!! limiting at cells adjacent to the ice front.
pure real function minmod2(a, b)
  real, intent(in) :: a !< First argument [Z ~> m]
  real, intent(in) :: b !< Second argument [Z ~> m]

  if (a > 0.0 .and. b > 0.0) then
    minmod2 = min(a, b)
  elseif (a < 0.0 .and. b < 0.0) then
    minmod2 = max(a, b)
  else
    minmod2 = 0.0
  endif
end function minmod2

!> Compute driving stress at B-grid nodes using a Pure DG(1) formulation.
!! Allows sub-element driving stress around grounding line
!! Evaluates the FEM weak-form integral using integration by parts.
!! To prevent B-grid null-space checkerboarding, the boundary integrals
!! utilize a Lax-Friedrichs (Rusanov) numerical flux. This ensures a
!! single-valued pressure at cell interfaces and applies a rigorous jump
!! penalty to aggressively damp grid-scale slope oscillations.
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
  real :: P_star           ! Lax-Friedrichs numerical flux at face [R Z L2 T-2 ~> kg s-2]
  real :: alpha_pen        ! Rusanov jump penalty coefficient [R L2 T-2 ~> kg m-1 s-2]
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
  real :: dsdx_gp, dsdy_gp  ! Surface slope at a qp [Z L-1 ~> nondim]
  real :: bottom_force_x, bottom_force_y ! Bed/water bottom drag forces [R Z L2 T-2 ~> kg s-2]
  real :: dphi_dx_ref, dphi_dy_ref ! Basis function derivatives in reference coordinates [nondim]
  real :: dphi_dx, dphi_dy  ! Basis function derivatives in physical coordinates [L-1 ~> m-1]
  real :: p_term_vol        ! Integrated-by-parts volume pressure term [R Z L2 T-2 ~> kg s-2]
  real :: scale             ! Multiplier applied to dsdx,dsdy to enforce max_surface_slope [nondim]
  real :: slope_mag         ! |grad(s)| magnitude at a qp [Z L-1 ~> nondim]
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
  real, dimension(2,2) :: slope_x_gp, slope_y_gp ! Per-QP surface slopes [Z L-1 ~> nondim]
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
            ISS%h_shelf(i,j), CS%h_x(i,j), CS%h_y(i,j), bed_corners, &
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

          h_gp = max(ISS%h_shelf(i,j) + ((CS%h_x(i,j)*xi_gp) + (CS%h_y(i,j)*eta_gp)), CS%min_h_shelf)
          dhdx_gp = CS%h_x(i,j) / a_qp
          dhdy_gp = CS%h_y(i,j) / d_qp

          bed_gp = ((bed_corners(1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                    (bed_corners(2,2) * (xquad(iq)   * xquad(jq))))  + &
                   ((bed_corners(2,1) * (xquad(iq)   * xquad(3-jq))) + &
                    (bed_corners(1,2) * (xquad(3-iq) * xquad(jq))))

          ! Reference-coord bed gradients
          dbdx_ref = (((bed_corners(1,1) * (-xquad(3-jq))) + &
                       (bed_corners(2,2) * ( xquad(jq))))  + &
                      ((bed_corners(2,1) * ( xquad(3-jq))) + &
                       (bed_corners(1,2) * (-xquad(jq)))))
          dbdy_ref = (((bed_corners(1,1) * (-xquad(3-iq))) + &
                       (bed_corners(2,2) * ( xquad(iq))))  + &
                      ((bed_corners(2,1) * (-xquad(iq))) + &
                       (bed_corners(1,2) * ( xquad(3-iq)))))
          dbdx_gp = dbdx_ref / a_qp
          dbdy_gp = dbdy_ref / d_qp

          if (CS%GL_couple) then
            if (CS%ground_frac(i,j)>0) then
              dsdx_gp = (1.0 - rhoi_rhow) * dhdx_gp
              dsdy_gp = (1.0 - rhoi_rhow) * dhdy_gp
              bottom_force_x = (rho**2 / rhow) * grav * h_gp * dhdx_gp
              bottom_force_y = (rho**2 / rhow) * grav * h_gp * dhdy_gp
            else
              dsdx_gp = dhdx_gp - dbdx_gp
              dsdy_gp = dhdy_gp - dbdy_gp
              bottom_force_x = rho * grav * h_gp * dbdx_gp
              bottom_force_y = rho * grav * h_gp * dbdy_gp
            endif
          else
            if (rhoi_rhow * h_gp - bed_gp <= 0.0) then
              dsdx_gp = (1.0 - rhoi_rhow) * dhdx_gp
              dsdy_gp = (1.0 - rhoi_rhow) * dhdy_gp
              bottom_force_x = (rho**2 / rhow) * grav * h_gp * dhdx_gp
              bottom_force_y = (rho**2 / rhow) * grav * h_gp * dhdy_gp
            else
              dsdx_gp = dhdx_gp - dbdx_gp
              dsdy_gp = dhdy_gp - dbdy_gp
              bottom_force_x = rho * grav * h_gp * dbdx_gp
              bottom_force_y = rho * grav * h_gp * dbdy_gp
            endif
          endif

          scale = 1.0
          if (CS%max_surface_slope > 0.0) then
            slope_mag = sqrt((dsdx_gp*dsdx_gp) + (dsdy_gp*dsdy_gp))
            scale = CS%max_surface_slope / max(slope_mag, CS%max_surface_slope)
          endif

          ! For slope diagnostics
          if (calc_slope_diag) then
            slope_x_gp(iq,jq) = dsdx_gp*scale
            slope_y_gp(iq,jq) = dsdy_gp*scale
          endif

          ! Weak-form Volume Integration by Parts
          p_term_vol = 0.5 * rho * grav * h_gp**2

          do n=1,2 ; do m=1,2
            phi_val = (merge(xquad(iq), xquad(3-iq), m == 2)) * &
                      (merge(xquad(jq), xquad(3-jq), n == 2))
            dphi_dx_ref = (merge(1.0, -1.0, m == 2)) * &
                          (merge(xquad(jq), xquad(3-jq), n == 2))
            dphi_dy_ref = (merge(xquad(iq), xquad(3-iq), m == 2)) * &
                          (merge(1.0, -1.0, n == 2))
            dphi_dx = dphi_dx_ref / a_qp
            dphi_dy = dphi_dy_ref / d_qp

            qp_dx(iq,jq,m,n) = scale * (weight * dphi_dx * p_term_vol + weight * phi_val * bottom_force_x)
            qp_dy(iq,jq,m,n) = scale * (weight * dphi_dy * p_term_vol + weight * phi_val * bottom_force_y)
          enddo ; enddo

        enddo ; enddo

        do n=1,2 ; do m=1,2
          vol_dx(m,n) = (qp_dx(1,1,m,n) + qp_dx(2,2,m,n)) + (qp_dx(1,2,m,n) + qp_dx(2,1,m,n))
          vol_dy(m,n) = (qp_dy(1,1,m,n) + qp_dy(2,2,m,n)) + (qp_dy(1,2,m,n) + qp_dy(2,1,m,n))
        enddo ; enddo
      endif

      if (calc_slope_diag) then
        if (.not. (CS%GL_regularize .and. CS%ground_frac(i,j) > 0.0 .and. CS%ground_frac(i,j) < 1.0)) then
          CS%sx_shelf(i,j) = 0.25*((slope_x_gp(1,1)+slope_x_gp(2,2)) + (slope_x_gp(1,2)+slope_x_gp(2,1)))
          CS%sy_shelf(i,j) = 0.25*((slope_y_gp(1,1)+slope_y_gp(2,2)) + (slope_y_gp(1,2)+slope_y_gp(2,1)))
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
      h_loc_A = max(ISS%h_shelf(i,j) + (((-0.5)*CS%h_x(i,j)) + ((-0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      h_loc_B = max(ISS%h_shelf(i,j) + (((-0.5)*CS%h_x(i,j)) + (( 0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      b_loc_A = bed_corners(1,1) ; b_loc_B = bed_corners(1,2)

      h_ngh_A = max(ISS%h_shelf(i-1,j) + ((( 0.5)*CS%h_x(i-1,j)) + ((-0.5)*CS%h_y(i-1,j))), CS%min_h_shelf)
      h_ngh_B = max(ISS%h_shelf(i-1,j) + ((( 0.5)*CS%h_x(i-1,j)) + (( 0.5)*CS%h_y(i-1,j))), CS%min_h_shelf)
      b_ngh_A = bed_corners(1,1) ; b_ngh_B = bed_corners(1,2)

      is_ext_bdry = ((CS%u_face_mask_bdry(I-1,j) == 2) .or. &
                    ((ISS%hmask(i-1,j) == 0 .or. ISS%hmask(i-1,j) == 2) .and. &
                     (CS%reentrant_x .or. (i+i_off /= gisc))))

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
          d_ocean = min(b_loc, rhoi_rhow * h_loc)
          P_star = 0.5 * grav * rhow * d_ocean**2
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Lax-Friedrichs Flux: {P} - 0.5*alpha*(h_right - h_left)
          ! Left is ngh, Right is loc
          alpha_pen = rho * grav * max(h_loc, h_ngh)
          P_star = 0.5 * (P_ngh + P_loc) - 0.5 * alpha_pen * (h_loc - h_ngh)
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
      h_loc_A = max(ISS%h_shelf(i,j) + ((( 0.5)*CS%h_x(i,j)) + ((-0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      h_loc_B = max(ISS%h_shelf(i,j) + ((( 0.5)*CS%h_x(i,j)) + (( 0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      b_loc_A = bed_corners(2,1) ; b_loc_B = bed_corners(2,2)

      h_ngh_A = max(ISS%h_shelf(i+1,j) + (((-0.5)*CS%h_x(i+1,j)) + ((-0.5)*CS%h_y(i+1,j))), CS%min_h_shelf)
      h_ngh_B = max(ISS%h_shelf(i+1,j) + (((-0.5)*CS%h_x(i+1,j)) + (( 0.5)*CS%h_y(i+1,j))), CS%min_h_shelf)
      b_ngh_A = bed_corners(2,1) ; b_ngh_B = bed_corners(2,2)

      is_ext_bdry = ((CS%u_face_mask_bdry(I,j) == 2) .or. &
                    ((ISS%hmask(i+1,j) == 0 .or. ISS%hmask(i+1,j) == 2) .and. &
                     (CS%reentrant_x .or. (i+i_off /= giec))))

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
          d_ocean = min(b_loc, rhoi_rhow * h_loc)
          P_star = 0.5 * grav * rhow * d_ocean**2
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Lax-Friedrichs Flux: {P} - 0.5*alpha*(h_right - h_left)
          ! Left is loc, Right is ngh
          alpha_pen = rho * grav * max(h_loc, h_ngh)
          P_star = 0.5 * (P_loc + P_ngh) - 0.5 * alpha_pen * (h_ngh - h_loc)
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
      h_loc_A = max(ISS%h_shelf(i,j) + (((-0.5)*CS%h_x(i,j)) + ((-0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      h_loc_B = max(ISS%h_shelf(i,j) + ((( 0.5)*CS%h_x(i,j)) + ((-0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      b_loc_A = bed_corners(1,1) ; b_loc_B = bed_corners(2,1)

      h_ngh_A = max(ISS%h_shelf(i,j-1) + (((-0.5)*CS%h_x(i,j-1)) + (( 0.5)*CS%h_y(i,j-1))), CS%min_h_shelf)
      h_ngh_B = max(ISS%h_shelf(i,j-1) + ((( 0.5)*CS%h_x(i,j-1)) + (( 0.5)*CS%h_y(i,j-1))), CS%min_h_shelf)
      b_ngh_A = bed_corners(1,1) ; b_ngh_B = bed_corners(2,1)

      is_ext_bdry = ((CS%v_face_mask_bdry(i,J-1) == 2) .or. &
                    ((ISS%hmask(i,j-1) == 0 .or. ISS%hmask(i,j-1) == 2) .and. &
                     (CS%reentrant_y .or. (j+j_off /= gjsc))))

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
          d_ocean = min(b_loc, rhoi_rhow * h_loc)
          P_star = 0.5 * grav * rhow * d_ocean**2
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Lax-Friedrichs Flux
          ! Left is ngh, Right is loc
          alpha_pen = rho * grav * max(h_loc, h_ngh)
          P_star = 0.5 * (P_ngh + P_loc) - 0.5 * alpha_pen * (h_loc - h_ngh)
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
      h_loc_A = max(ISS%h_shelf(i,j) + (((-0.5)*CS%h_x(i,j)) + (( 0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      h_loc_B = max(ISS%h_shelf(i,j) + ((( 0.5)*CS%h_x(i,j)) + (( 0.5)*CS%h_y(i,j))), CS%min_h_shelf)
      b_loc_A = bed_corners(1,2) ; b_loc_B = bed_corners(2,2)

      h_ngh_A = max(ISS%h_shelf(i,j+1) + (((-0.5)*CS%h_x(i,j+1)) + ((-0.5)*CS%h_y(i,j+1))), CS%min_h_shelf)
      h_ngh_B = max(ISS%h_shelf(i,j+1) + ((( 0.5)*CS%h_x(i,j+1)) + ((-0.5)*CS%h_y(i,j+1))), CS%min_h_shelf)
      b_ngh_A = bed_corners(1,2) ; b_ngh_B = bed_corners(2,2)

      is_ext_bdry = ((CS%v_face_mask_bdry(i,J) == 2) .or. &
                    ((ISS%hmask(i,j+1) == 0 .or. ISS%hmask(i,j+1) == 2) .and. &
                     (CS%reentrant_y .or. (j+j_off /= gjec))))

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
          d_ocean = min(b_loc, rhoi_rhow * h_loc)
          P_star = 0.5 * grav * rhow * d_ocean**2
        else
          h_ngh = (1.0 - t_face)*h_ngh_A + t_face*h_ngh_B
          b_ngh = (1.0 - t_face)*b_ngh_A + t_face*b_ngh_B
          if (rhoi_rhow * h_ngh - b_ngh > 0.0) then
            P_ngh = 0.5 * grav * rho * h_ngh**2
          else
            P_ngh = 0.5 * grav * (1.0 - rhoi_rhow) * rho * h_ngh**2
          endif

          ! Lax-Friedrichs Flux
          ! Left is loc, Right is ngh
          alpha_pen = rho * grav * max(h_loc, h_ngh)
          P_star = 0.5 * (P_loc + P_ngh) - 0.5 * alpha_pen * (h_ngh - h_loc)
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

!> Subgrid GL-band volume integral of the driving stress for the DG path.
!! Evaluates the unified integration-by-parts weak form over nsub x nsub sub-cells.
subroutine calc_shelf_driving_stress_DG_subgrid(CS, Phisub, &
    h_shelf_cell, h_x_cell, h_y_cell, bed_corners, &
    dxCv_S, dxCv_N, dyCu_W, dyCu_E, &
    rho, rhow, rhoi_rhow, grav, vol_dx, vol_dy, sx_shelf, sy_shelf, calc_slope_diag)
  type(ice_shelf_dyn_CS), intent(in) :: CS    !< Ice shelf control structure
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Sub-grid quadrature weights [nondim]
  real, intent(in) :: h_shelf_cell   !< Cell-averaged ice thickness [Z ~> m]
  real, intent(in) :: h_x_cell       !< DG x-slope moment [Z ~> m]
  real, intent(in) :: h_y_cell       !< DG y-slope moment [Z ~> m]
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
  logical :: calc_slope_diag ! True if slope diagnostics will be calculated

  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: contr_sub_dx, contr_sub_dy
  real, dimension(SIZE(Phisub,3),SIZE(Phisub,3),2,2) :: slope_x_gp, slope_y_gp ! Per-QP surf slopes [Z L-1 ~> nondim]
  real, dimension(2,2) :: slope_x, slope_y ! slope sums for qps with same position within subcells [Z L-1 ~> nondim]
  real, dimension(2,2,2,2) :: qp_dx, qp_dy
  real :: xi_sub, eta_sub   ! DG reference coords at sub-qp ([-0.5,0.5]) [nondim]
  real :: h_gp              ! Ice thickness at sub-qp [Z ~> m]
  real :: dhdx_gp, dhdy_gp  ! Thickness gradients at sub-qp, physical coords [Z L-1 ~> nondim]
  real :: bed_gp            ! Bed elevation at sub-qp [Z ~> m]
  real :: dbdx_gp, dbdy_gp  ! Bed gradients at sub-qp, physical coords [Z L-1 ~> nondim]
  real :: dbdx_ref, dbdy_ref ! Bed gradients in reference coords (pre-Jacobian) [Z ~> m]
  real :: dsdx_gp, dsdy_gp   ! Surface slope at sub-qp [Z L-1 ~> nondim]
  real :: bottom_force_x, bottom_force_y ! Bed/water bottom drag forces [R Z L2 T-2 ~> kg s-2]
  real :: dphi_dx_ref, dphi_dy_ref ! Basis function derivatives in reference coordinates [nondim]
  real :: dphi_dx, dphi_dy  ! Basis function derivatives in physical coordinates [L-1 ~> m-1]
  real :: y_marginal_1, y_marginal_2, x_marginal_1, x_marginal_2 ! Marginal sums [nondim]
  real :: p_term_vol        ! Integrated-by-parts volume pressure term [R Z L2 T-2 ~> kg s-2]
  real :: a, d              ! Per-sub-qp interpolated cell-edge spacings [L ~> m]
  real :: weight            ! Per-sub-qp quadrature weight [L2 ~> m2]
  real :: subarea           ! 1/nsub^2 [nondim]
  real :: scale, slope_mag  ! max_surface_slope clamp [nondim]
  integer :: nsub, i, j, qx, qy, m, n

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  do j=1,nsub ; do i=1,nsub
    qp_dx(:,:,:,:) = 0.0 ; qp_dy(:,:,:,:) = 0.0
    do qy=1,2 ; do qx=1,2
      ! Reference coords at the sub-qp: xi_sub = a_right(qx,i) - 0.5; marginal
      ! sum of Phisub over k for l=2 gives a_right(qx,i).
      xi_sub  = (Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)) - 0.5
      eta_sub = (Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)) - 0.5

      h_gp = max(h_shelf_cell + ((h_x_cell*xi_sub) + (h_y_cell*eta_sub)), CS%min_h_shelf)

      ! Bed at sub-qp: rotation-paired Phisub contraction of bed_corners.
      bed_gp = ((Phisub(qx,qy,i,j,1,1)*bed_corners(1,1)) + (Phisub(qx,qy,i,j,2,2)*bed_corners(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*bed_corners(1,2)) + (Phisub(qx,qy,i,j,2,1)*bed_corners(2,1)))

      ! Per-sub-qp metric via Phisub marginal sums (same pattern as CG_action_subgrid_basal).
      a = (dxCv_S * (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1))) + &
          (dxCv_N * (Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)))
      d = (dyCu_W * (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2))) + &
          (dyCu_E * (Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)))
      weight = 0.25 * subarea * (a * d)

      y_marginal_1 = Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1)
      y_marginal_2 = Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2)
      x_marginal_1 = Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2)
      x_marginal_2 = Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)

      dhdx_gp = h_x_cell / a
      dhdy_gp = h_y_cell / d

      ! Reference-coord bed gradients: derivative of bilinear corner-basis at sub-qp.
      dbdx_ref = (((bed_corners(1,1) * (-(Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1)))) + &
                   (bed_corners(2,2) * ( (Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2))))) + &
                  ((bed_corners(2,1) * ( (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,2,1)))) + &
                   (bed_corners(1,2) * (-(Phisub(qx,qy,i,j,1,2) + Phisub(qx,qy,i,j,2,2))))))
      dbdy_ref = (((bed_corners(1,1) * (-(Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2)))) + &
                   (bed_corners(2,2) * ( (Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2))))) + &
                  ((bed_corners(2,1) * (-(Phisub(qx,qy,i,j,2,1) + Phisub(qx,qy,i,j,2,2)))) + &
                   (bed_corners(1,2) * ( (Phisub(qx,qy,i,j,1,1) + Phisub(qx,qy,i,j,1,2))))))
      dbdx_gp = dbdx_ref / a
      dbdy_gp = dbdy_ref / d

      if (rhoi_rhow * h_gp - bed_gp <= 0.0) then
        ! Floating: bottom force is water pressure on the sloped draft
        dsdx_gp = (1.0 - rhoi_rhow) * dhdx_gp
        dsdy_gp = (1.0 - rhoi_rhow) * dhdy_gp
        bottom_force_x = (rho**2 / rhow) * grav * h_gp * dhdx_gp
        bottom_force_y = (rho**2 / rhow) * grav * h_gp * dhdy_gp
      else
        ! Grounded: bottom force is bed pressure on the sloped bed
        bottom_force_x = rho * grav * h_gp * dbdx_gp
        bottom_force_y = rho * grav * h_gp * dbdy_gp
        dsdx_gp = dhdx_gp - dbdx_gp
        dsdy_gp = dhdy_gp - dbdy_gp
      endif

      if (CS%max_surface_slope > 0.0) then
        slope_mag = sqrt((dsdx_gp*dsdx_gp) + (dsdy_gp*dsdy_gp))
        scale = CS%max_surface_slope / max(slope_mag, CS%max_surface_slope)
      endif

      ! For slope diagnostics
      if (calc_slope_diag) then
        slope_x_gp(i,j,qx,qy) = dsdx_gp*scale
        slope_y_gp(i,j,qx,qy) = dsdy_gp*scale
      endif

      ! Unified Weak-form Volume Integration applied to subgrid
      p_term_vol = 0.5 * rho * grav * h_gp**2

      do n=1,2 ; do m=1,2
        dphi_dx_ref = merge(1.0, -1.0, m==2) * merge(y_marginal_2, y_marginal_1, n==2)
        dphi_dy_ref = merge(x_marginal_2, x_marginal_1, m==2) * merge(1.0, -1.0, n==2)
        dphi_dx = dphi_dx_ref / a
        dphi_dy = dphi_dy_ref / d

        qp_dx(qx,qy,m,n) = scale * (weight * dphi_dx * p_term_vol + weight * Phisub(qx,qy,i,j,m,n) * bottom_force_x)
        qp_dy(qx,qy,m,n) = scale * (weight * dphi_dy * p_term_vol + weight * Phisub(qx,qy,i,j,m,n) * bottom_force_y)
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
    enddo; enddo

    sx_shelf = 0.25*((slope_x(1,1)+slope_x(2,2)) + (slope_x(1,2)+slope_x(2,1)))/nsub
    sy_shelf = 0.25*((slope_y(1,1)+slope_y(2,2)) + (slope_y(1,2)+slope_y(2,1)))/nsub
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

end module MOM_ice_shelf_dynamics
