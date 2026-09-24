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
use MOM_domains, only : FOLD_NORTH_EDGE
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
use MOM_ice_shelf_initialize, only : corner_cell_weights, nodal_cell_mean, cell_face_lengths
implicit none ; private

#include <MOM_memory.h>

public register_ice_shelf_dyn_restarts, initialize_ice_shelf_dyn, update_ice_shelf, IS_dynamics_post_data
public ice_time_step_CFL, ice_shelf_dyn_end, change_in_draft, write_ice_shelf_energy
public shelf_advance_front, ice_shelf_min_thickness_calve, calve_to_mask, volume_above_floatation
public reset_DG_to_cellmean_at_cell, zero_DG_h_nodal
public DG_nodal_thickness_ptr
public accumulate_DG_source_rate
public calc_prescribed_basal_melt
public masked_var_grounded

! SSA inner solver flags
integer, parameter :: INNER_CG = 1       !< Conjugate gradient (default)
integer, parameter :: INNER_MINRES = 2   !< MINRES
integer, parameter :: INNER_CR = 3       !< Conjugate residual

! DG(1) mode-damper gate forms, selected by DG1_TILT_DAMP_GATE
integer, parameter :: DAMP_GATE_THICKNESS = 0 !< Gate the zigzag against the ice thickness
integer, parameter :: DAMP_GATE_SLOPE = 1     !< Gate it against the real surface slope
integer, parameter :: DAMP_GATE_AGREEMENT = 2 !< Gate it on how well the readings agree

! DG(1) mode-damper reductions, selected by DG1_TILT_DAMP_REDUCE
integer, parameter :: DAMP_RED_MINMOD = 0 !< Compare every reading with every other one
integer, parameter :: DAMP_RED_PAIRED = 1 !< Add the two one-sided readings before comparing
                                          !! them with the centred one

! DG(1) mode-damper detectors, selected by DG1_TILT_DAMP_DETECTOR
integer, parameter :: DAMP_EDGE_OFF = 0 !< Leave a cell that has lost a neighbour to the
                                        !! stencil readings, which cannot cancel there
integer, parameter :: DAMP_EDGE_REF = 1 !< Take the estimates from the cell means instead, and
                                        !! stand down where the run crosses a flotation break
integer, parameter :: DAMP_EDGE_REF_GL = 2 !< The same, but read the means across a flotation
                                        !! break as well

! Which fields the agreement detector reads, selected by DG1_TILT_DAMP_FIELDS
integer, parameter :: DAMP_FLD_BOTH = 0 !< The thickness and the surface form must agree
integer, parameter :: DAMP_FLD_SURF = 1 !< The surface form alone
integer, parameter :: DAMP_FLD_THCK = 2 !< The thickness alone
integer, parameter :: DAMP_FLD_GL_TH = 3 !< Both, except at a cell the flotation contour
                                         !! crosses, where the surface form is a fraction of
                                         !! the bed removed and means nothing
integer, parameter :: DAMP_DET_TODAY = 0 !< Tilt Laplacian against bed and mean-supported floors
integer, parameter :: DAMP_DET_AGREE = 1 !< Stencil agreement over three stencils and two fields
integer, parameter :: DAMP_DET_HYBRID = 2 !< Each cell chooses: the tilt Laplacian where the
                                          !! stencil stays on one side of the flotation
                                          !! contour, stencil agreement where it does not


! Sentinel returned by the Coulomb fB routines when the effective pressure is zero. The sliding law
! tau_b = C |u|^m / (1 + fB |u|^q)^m gives exactly zero drag as N -> 0, but it encodes that limit as
! fB -> infinity, which is not representable. Any negative fB therefore means "no Coulomb drag here",
! and compute_basal_coef takes that branch instead of evaluating the divergent expression. This is
! what lets CF_MinN be set to zero: with a positive CF_MinN the effective pressure is floored before
! fB is formed and the sentinel is never produced.
real, parameter :: FB_NO_COULOMB_DRAG = -1.0 !< fB value meaning zero effective pressure [(T L-1)^CF_PostPeak]

!> Flotation-branch amplification in the DG(1) artificial-viscosity jump-mode rate
!! lambda = 4*amp*c*u_eff/dx_perp. The harmonic slope mean of dg1_wb_slope_mean
!! makes amp exactly 2 at every face [nondim].
real, parameter :: DG1_WB_JUMP_RATE_AMP = 2.0

! CISM-style grounding-line treatment modes (CISM_FRICTION, CISM_TAUD)
integer, parameter :: CISM_OFF = 0       !< No CISM-style treatment of this term
integer, parameter :: CISM_LOCAL = 1     !< Quadrant grounded fraction with the local/lumped nodal assembly
integer, parameter :: CISM_INTEGRATE = 2 !< Quadrant grounded fraction with the consistent element assembly

! SEP2 sub-element quadrature constants (GROUNDING_LINE_SUBGRID_SCHEME="SEP2").
real, parameter :: SEP2_W23 = 2.0/3.0    !< Heavy vertex weight of the interior 3-pt triangle rule [nondim]
real, parameter :: SEP2_W16 = 1.0/6.0    !< Light vertex weight of the interior 3-pt triangle rule [nondim]
real, parameter :: SEP2_TRI3 = 0.25/3.0  !< Per-QP reference measure of a whole parent triangle [nondim]
real, parameter, dimension(2) :: SEP2_GP = (/ 0.5*(1.0 - sqrt(1.0/3.0)), &
                                              0.5*(1.0 + sqrt(1.0/3.0)) /)
                                         !< 2-pt Gauss abscissae on [0,1] [nondim].
                                         !! Built from this expression rather than from decimal
                                         !! literals because the pair must satisfy both
                                         !! 1-GP(1) == GP(2) and 1-GP(2) == GP(1) to the last bit.
                                         !! A quarter turn of the grid maps each abscissa to the
                                         !! complement of the other, so a pair that is asymmetric
                                         !! by one ulp makes every quadrature built on it
                                         !! irreproducible under rotation.
real, parameter, dimension(2) :: SEP2_GC = (/ SEP2_GP(2), SEP2_GP(1) /)
                                         !< Complementary Gauss factors (1-abscissa), stored as the
                                         !! same values swapped so reflection orbits are exact [nondim]

! TVD slope limiters for thickness advection (ICE_SHELF_ADVECT_LIMITER)
integer, parameter :: LIMITER_VANLEER = 0   !< Van Leer limiter (original scheme)
integer, parameter :: LIMITER_SUPERBEE = 1  !< Superbee limiter (least diffusive; STREAMICE default)
integer, parameter :: LIMITER_MINMOD = 2    !< Minmod limiter (most diffusive)
integer, parameter :: LIMITER_MC = 3        !< Monotonized-central limiter (between Van Leer and superbee)

! How the DG(1) nodal thickness source is shared between cells meeting at a corner.
! All are exactly mass-conservative.
integer, parameter :: SRC_OP_AVERAGED = 0 !< Corner takes the cell_mean_w-weighted average of its cells
integer, parameter :: SRC_OP_LOCAL = 1  !< Corner takes its own cell's rate; nothing crosses a face
integer, parameter :: SRC_OP_SUBGRID = 2 !< Averaged over the floating part of each corner's support
                                        !! and returned by floating fraction. Needs xi_basal.

! Grounding-line treatment of the prescribed ice-only basal melt (ICE_ONLY_BASAL_MELT_GLP).
! Named after Leguy, Lipscomb & Asay-Davis (2021) sec. 2.3 and Seroussi & Morlighem (2018) sec. 2.
integer, parameter :: MELT_GLP_FMP = 0  !< Full melt: the fully-floating rate is applied in every
                                        !! ice-covered cell regardless of its grounded fraction.
integer, parameter :: MELT_GLP_FCMP = 1 !< Flotation-condition melt: full rate where the cell centre
                                        !! satisfies the flotation condition, zero otherwise.
integer, parameter :: MELT_GLP_PMP = 2  !< Partial melt: the rate is scaled by the floating area
                                        !! fraction of the cell. Equivalent to Seroussi's SEM1.
integer, parameter :: MELT_GLP_NMP = 3  !< No melt: zero rate in every partly grounded cell.
integer, parameter :: MELT_GLP_SEM2 = 4 !< Seroussi & Morlighem (2018) SEM2: PMP's cell total,
                                        !! distributed by nodal floating fraction. DG only.


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
  real, pointer, dimension(:,:,:,:) :: h_nodal => NULL() !< DG(1) Q1 thickness at each cell's own corners,
                                                       !! h_nodal(i,j,a,b) with a = W/E, b = S/N [Z ~> m].
  real, pointer, dimension(:,:,:,:,:,:) :: Minv_nodal => NULL() !< Inverse Q1 mass matrix of cell (i,j),
                                   !! Minv_nodal(i,j,a,b,p,q) = Minv_xi(a,p)*Minv_eta(b,q). It maps the
                                   !! right-hand side at corner (p,q) to the rate at corner (a,b);
                                   !! a and p are 1 (W) or 2 (E), b and q are 1 (S) or 2 (N) [L-2 ~> m-2].
  real, pointer, dimension(:,:,:,:) :: cell_mean_w => NULL() !< Integral over the cell of the basis function of
                                   !! corner (a,b). The four sum to the cell area, and sum w*h / sum w is
                                   !! the area-weighted cell mean of a corner field. Used for cell means,
                                   !! nodal source projection and mass accounting [L2 ~> m2].
  real, pointer, dimension(:,:) :: h_source_rate => NULL() !< Cell-mean thickness source (basal + surface)
                                                       !! accumulated since the last DG advect step [Z T-1 ~> m s-1].
  real, pointer, dimension(:,:) :: h_source_rate_bmb => NULL() !< Basal part of h_source_rate [Z T-1 ~> m s-1].
  real, pointer, dimension(:,:) :: h_source_rate_last => NULL() !< h_source_rate used by the last DG advect
                                                       !! step, kept for diagnostics [Z T-1 ~> m s-1].
  real, pointer, dimension(:,:) :: C_basal_friction => NULL()!< Coefficient in sliding law tau_b = C u^(n_basal_fric),
                               !! units of [R L Z T-2 (s m-1)^(n_basal_fric) ~> Pa (s m-1)^(n_basal_fric)]
  real, pointer, dimension(:,:) :: coef_prefactor => NULL() !< Pre-computed area*C_basal_friction*L_T_to_m_s for
                               !! basal friction quadrature evaluation [R L2 Z T-1 ~> kg s-1].
  real, pointer, dimension(:,:) :: fB_elem => NULL()        !< Pre-computed element-level Coulomb fB parameter
                               !! [(T L-1)^CF_PostPeak]; 0 for Weertman.
                               !! Updated each outer iteration by calc_shelf_basal_prefactors.
  real, pointer, dimension(:,:) :: coef_prefactor_node => NULL() !< Pre-computed area_node*C_node*L_T_to_m_s at
                               !! B-grid nodes for the local (CISM_FRICTION='local') diagonal drag,
                               !! C_node an area-weighted 4-cell average and area_node the ice-restricted
                               !! nodal control volume [R L2 Z T-1 ~> kg s-1].
  real, pointer, dimension(:,:) :: fB_node => NULL()        !< Pre-computed nodal Coulomb fB parameter at B-grid
                               !! nodes for CISM_FRICTION='local' [(T L-1)^CF_PostPeak]; 0 for Weertman.
  real, pointer, dimension(:,:) :: area_node => NULL()      !< Nodal control-volume area for the local
                               !! (CISM_FRICTION='local') drag: the sum of the surrounding cells'
                               !! lumped corner areas (0.25*areaT each), i.e. the same control volume the
                               !! lumped driving stress and the CG_action element assembly integrate over.
                               !! With LOCAL_NODE_FULL_AREA the sum is over all four in-domain cells (the
                               !! CISM convention, dx*dy at every active vertex); otherwise it is restricted
                               !! to the ice-covered cells. Either way it halves at domain-edge nodes, so the
                               !! local taud and friction share one control volume and the wall-node x-force
                               !! balance stays meridionally symmetric [L2 ~> m2].
  real :: alpha_coulomb = 1.0  !< Coulomb prefactor (CF_PostPeak-1)^(CF_PostPeak-1)/CF_PostPeak^CF_PostPeak [nondim]
  real, pointer, dimension(:,:) :: OD_rt => NULL()         !< A running total for calculating OD_av [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac_rt => NULL() !< A running total for calculating ground_frac.
  real, pointer, dimension(:,:) :: OD_av => NULL()         !< The time average open ocean depth [Z ~> m].
  real, pointer, dimension(:,:) :: ground_frac => NULL()   !< Fraction of the time a cell is "exposed", i.e. the column
                               !! thickness is below a threshold and interacting with the rock [nondim].  When this
                               !! is 1, the ice-shelf is grounded
  real, pointer, dimension(:,:) :: f_ground_node => NULL() !< Analytic grounded ice fraction at B-grid
                               !! nodes (vertices) from the quadrant grounding-line parameterization
                               !! (Leguy et al. 2021). Multiplies basal friction when CISM_FRICTION
                               !! is set. 1 = fully grounded, 0 = fully floating [nondim].
  real, pointer, dimension(:,:) :: f_ground_cell => NULL() !< Analytic grounded ice fraction at cell
                               !! centers from the same quadrant parameterization (shares the per-cell
                               !! quadrant areas with f_ground_node, so the two grids carry mutually
                               !! consistent grounded areas). Used to blend the surface for the FV
                               !! driving stress [nondim].
  real, pointer, dimension(:,:) :: H_node => NULL() !< The ice shelf thickness at B-grid corners,
                               !! set by interpolate_H_to_B in update_grounded_geometry and used by
                               !! the sub-element basal friction in the velocity solve.  It is only
                               !! nonzero with GROUNDING_LINE_INTERPOLATE and without DG thickness,
                               !! which are the cases that read it [Z ~> m].
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
  logical :: grounded_geom_current = .false. !< If true, the grounded geometry fields set by
                               !! update_grounded_geometry have already been refreshed for the current
                               !! ice thickness, so the velocity solve does not have to repeat the work.

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
  logical :: FV_GL_one_sided !< If true, use one-sided finite-volume differences to evaluate the
                            !! driving stress in the cells on either side of the grounding line,
                            !! following Cornford et al. (2013) eqs 27-29, rather than the
                            !! centered difference (their eq 25) that would straddle the
                            !! grounding line. Only used by the FV (non-DG) driving stress.

  logical :: ice_only_basal_melt !< If true, the ice-only (solo) driver applies a prescribed
                            !! depth-dependent basal melt rate under floating ice, following
                            !! Leguy et al. (2021) eq. 18 (= Seroussi & Morlighem 2018 eq. 4).
                            !! Has no effect in coupled runs, where the melt rate comes from the
                            !! ocean via shelf_calc_flux.
  integer :: ice_only_melt_glp !< The grounding-line treatment of the prescribed ice-only basal
                            !! melt, one of MELT_GLP_FMP, MELT_GLP_FCMP, MELT_GLP_PMP or
                            !! MELT_GLP_NMP.
  real :: ice_only_melt_scale !< A factor multiplying the whole prescribed ice-only basal melt
                            !! profile, so 5.0 gives the high-melt experiments of Leguy et al.
                            !! (2021) sec. 4.3 [nondim].

  logical :: gl_quad_friction !< If true, scale basal friction by an analytic nodal grounded
                            !! fraction (CS%f_ground_node) computed by the quadrant grounding-line
                            !! parameterization of Leguy et al. (2021), replacing the geometric
                            !! sub-cell (Phisub) friction integration in grounding-line cells.
  logical :: use_sep2       !< If true (GROUNDING_LINE_SUBGRID_SCHEME="SEP2"), grounding-line cells
                            !! are split geometrically into grounded/floating sub-elements
                            !! (Seroussi et al. 2014 SEP2, extended to quadrilaterals) for the basal
                            !! friction, driving stress, and grounded fraction, instead of the
                            !! uniform Phisub sub-sampling ("SEP3").
  logical :: fv_taud_vertex_grad !< If true, the FV (non-DG) driving stress evaluates the surface
                            !! gradient directly at B-grid nodes from the four surrounding cell
                            !! centers (Lipscomb et al. 2019 eq. 14, "option 3" margins), instead of
                            !! the wider cell-centroid centered difference. Less smeared across the
                            !! grounding line. Mutually exclusive with FV_GL_ONE_SIDED_TAUD.
  integer :: cism_friction = CISM_OFF !< CISM-style grounding-line treatment of the basal friction,
                            !! set by CISM_FRICTION; drives gl_quad_friction and local_basal_friction.
  integer :: cism_taud = CISM_OFF !< CISM-style grounding-line treatment of the FV driving stress,
                            !! set by CISM_TAUD; drives fv_taud_vertex_grad and local_fv_taud_vertex.
  logical :: local_fv_taud_vertex !< If true (CISM_TAUD='local'), assemble the
                            !! nodal driving stress by the local/lumped method (Lipscomb 2019 A4; CISM
                            !! HO_ASSEMBLE_TAUD_LOCAL): tau_d at a node uses that node's slope alone over
                            !! its nodal control mass. If false, use the consistent element-quadrature
                            !! assembly. Local is the CISM-faithful default (co-locates with a local
                            !! basal friction); consistent matches a consistent-mass friction.
  logical :: local_basal_friction !< If true, assemble basal drag with a local/nodal diagonal (CISM
                            !! HO_ASSEMBLE_BETA_LOCAL): drag at a node = beta(node)*areaBu*u(node),
                            !! beta from a nodal C, nodal velocity, and the nodal grounded fraction
                            !! f_ground_node, with no element integration or neighbor coupling. Requires
                            !! CISM_FRICTION (for f_ground_node). Pairs with CISM_TAUD='local'
                            !! to reproduce the all-local CISM/Leguy-2021 grounding-line setup.
  logical :: i2n_friction !< If true, the FV (non-DG) basal friction and the Coulomb
                            !! effective pressure are integrated over the sub-element grounding-line
                            !! partition selected by GROUNDING_LINE_SUBGRID_SCHEME, using the corner
                            !! thickness and flotation fields CS%H_corner and CS%fls_corner instead of
                            !! corner H with a cell-constant bed. Requires GROUNDING_LINE_INTERPOLATE.
  logical :: i2n_taud !< If true, the FV (non-DG) driving stress is integrated over the same
                            !! sub-element partition used by I2N_FRICTION, with the surface
                            !! elevation reconstructed as S = (1-r)*H + max(fls,0) so that its slope
                            !! kink lies exactly on the partition's grounding line. Requires
                            !! GROUNDING_LINE_INTERPOLATE.
  real, pointer, dimension(:,:) :: H_corner => NULL() !< Ice thickness interpolated to B-grid corners
                            !! with dual-cell Lagrange weights over the included cells only
                            !! (FV_SUBGRID_GL_* paths) [Z ~> m].
  real, pointer, dimension(:,:) :: fls_corner => NULL() !< Flotation deficit r*h - bed interpolated to
                            !! B-grid corners with the same weights and cell set as CS%H_corner, so
                            !! that fls = r*H - bed holds pointwise at every corner [Z ~> m].
  logical, pointer, dimension(:,:) :: corner_valid => NULL() !< True where CS%H_corner and
                            !! CS%fls_corner have at least one contributing cell.
  real, pointer, dimension(:,:,:) :: corner_wt => NULL() !< Dual-cell Q1 (Lagrange) interpolation
                            !! weights of the 4 cells surrounding each B-grid node, ordered
                            !! SW, SE, NW, NE; each cell is weighted by the opposite cell's spacing,
                            !! so a linear field is reproduced exactly at the node. Equal to 1/4 on a
                            !! uniform Cartesian grid [nondim].
  integer :: adv_thickness_limiter = LIMITER_VANLEER !< TVD slope limiter used for thickness
                            !! advection in ice_shelf_advect_thickness_x/y (LIMITER_VANLEER,
                            !! LIMITER_SUPERBEE, LIMITER_MINMOD, or LIMITER_MC).
  logical :: adv_cfl_weight = .false. !< If true, weight the thickness-advection slope reconstruction by the
                            !! Lax-Wendroff (1-CFL) factor (as in STREAMICE), for a time-accurate
                            !! 2nd-order flux. If false, use the full spatial slope (original behavior).

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
  logical :: use_DG_thickness     !< If true, use the DG(1) nodal thickness with unsplit advection.
  logical :: init_bed_nodal       !< If true, read the bed at B-grid nodes into CS%bed_node and set
                                  !! CS%bed_elev to its cell means. Required by USE_DG_THICKNESS.
  integer :: dg_basal_source_op   !< Sets the SRC_OP_* method (LOCAL, AVERAGED, or SUBGRID) for how much of a cell's
                                  !! basal melt is shared with neighbors it meets at a corner
  real, pointer, dimension(:,:,:,:) :: xi_basal => NULL() !< Floating share of corner (a,b)'s support
                                  !! from the active sub-element partition: 1 fully floating, 0 fully
                                  !! grounded. 1 where no sub-element partition applies [nondim].
  logical :: dg_basal_source_sem2 !< If true, distribute the basal DG(1) source by CS%xi_basal (SEM2).
  logical :: dg_surface_source_local !< If true, apply the surface DG(1) source with SRC_OP_LOCAL,
                                  !! otherwise with SRC_OP_AVERAGED.
  ! DG(1) tilt/twist mode damper (dg_nodal_mode_damp_rate)
  logical :: dg_tilt_damp         !< If true, the mode damper damps the grid-scale in-cell tilt modes.
  real :: dg_tilt_damp_r_hi       !< Mode damper: normalized detector excess at which the gate is
                                  !! fully open [nondim].
  logical :: dg_twist_damp        !< If true, the mode damper damps the grid-scale in-cell twist mode.
  integer :: dg_damp_gate_form    !< Mode damper: which gate scales the damping, one of
                                  !! DAMP_GATE_THICKNESS, DAMP_GATE_SLOPE or DAMP_GATE_AGREEMENT.
  integer :: dg_damp_detector     !< Mode damper: which detector finds the zigzag, one of
                                  !! DAMP_DET_TODAY or DAMP_DET_AGREE.
  integer :: dg_damp_edge_rule   !< Mode damper: what to do at a cell that has lost a
                                 !! neighbour; see the DG1_TILT_DAMP_EDGE parameter.
  logical :: dg_damp_centred_amp !< Mode damper: if true, the minmod of the readings decides
                                 !! whether to damp and the centred reading sets how much.
  integer :: dg_damp_reduce      !< Mode damper: how the agreement detector reduces its
                                 !! readings to one value; see DG1_TILT_DAMP_REDUCE.
  integer :: dg_damp_fields      !< Mode damper: which fields the agreement detector reads;
                                 !! see DG1_TILT_DAMP_FIELDS.
  logical :: dg_damp_single_rule  !< Mode damper: if true, the agreement detector leaves a slope
                                  !! alone in a cell that has one stencil when that stencil
                                  !! crosses the grounding line.
  logical :: dg_damp_filter       !< Mode damper: if true, high-pass the removal so that only
                                  !! its grid-scale part is taken away.
  real :: dg_damp_rho_g           !< Mode damper: disagreement between the detector readings,
                                  !! relative to the part they agree on, at which the agreement
                                  !! gate is half open [nondim].
  real :: dg_damp_rho_s           !< Mode damper: zigzag surface slope, relative to the real surface
                                  !! slope, at which the slope gate is fully open [nondim].
  real :: dg_damp_slope_floor     !< Mode damper: smallest real surface slope the slope gate accepts
                                  !! as a reference [nondim].
  real, pointer, dimension(:,:) :: dg_damp_ahat => NULL() !< Mode damper: L2 cell norm of the
                                  !! detector output over the three modes [Z ~> m].
  real, pointer, dimension(:,:) :: dg_damp_spread => NULL() !< Mode damper: largest spread of the
                                  !! detector readings over the three modes [Z ~> m].
  real, pointer, dimension(:,:) :: dg_damp_tend => NULL() !< Mode damper: L2 cell norm of the
                                  !! correction rate [Z T-1 ~> m s-1].
  real, pointer, dimension(:,:) :: dg_damp_gate => NULL() !< Mode damper: largest gate over the
                                  !! three modes [nondim].
  real, pointer, dimension(:,:) :: dg_damp_want => NULL() !< Mode damper: rate requested, summed over
                                  !! modes, before the transport credit [T-1 ~> s-1].
  real, pointer, dimension(:,:) :: dg_damp_got => NULL() !< Mode damper: rate delivered, summed over
                                  !! modes [T-1 ~> s-1].
  logical :: dg_damp_advective    !< Mode damper: if true, subtract the removal rate the upwind
                                  !! transport already gives.
  real :: dg_damp_advective_c     !< Mode damper: coefficient c in the rates c*U_CUT/dx and
                                  !! c*|u_n|/dx [nondim].
  real :: dg_damp_u_cut           !< Mode damper: speed whose upwind removal rate c*u_cut/dx the
                                  !! damper supplies, less what the flow gives [L T-1 ~> m s-1].
  real :: dg_damp_kink_tol        !< Mode damper: grounded-fraction misfit at which the twist's kink
                                  !! reconstruction loses all confidence [nondim].
  real :: dg_damp_kink_fit_tol    !< Mode damper: relative slope misfit at which the tilt's kink
                                  !! reconstruction loses all confidence; <= 0 skips it [nondim].
  logical :: dg_damp_kink_ref     !< Mode damper: if true, reconstruct the slope break at the flotation
                                  !! contour instead of exempting the cell.
  logical :: dg_damp_excess_only  !< Mode damper: if true, remove only the detector part no reference
                                  !! explains.
  integer :: dg_damp_gl_reach     !< Mode damper grounding-line protection reach: 0 bisected cell,
                                  !! 1 detector stencil, 2 full reference stencil.
  logical :: dg_tilt_damp_dt_warned = .false. !< Mode damper: true once the short-relaxation warning
                                  !! is issued.
  ! DG(1) artificial viscosity (DG1_nodal_spatial_operator)
  real :: dg_art_visc_advect_coef !< Artificial viscosity: coefficient on |u_face| in u_eff [nondim].
  real :: dg_art_visc_strain_coef !< Artificial viscosity: coefficient on eps_e_face*dx_perp in u_eff [nondim].
  real :: dg_art_visc_advect_L_ref !< Artificial viscosity: if positive, scale the |u_face| term by
                                  !! dx_perp/L_ref so its decay rate is resolution-independent [L ~> m].
  real :: dg_art_visc_tau_floor   !< Artificial viscosity: if positive, add dx_perp/tau_floor to u_eff [T ~> s].
  logical :: dg_tilt_relax        !< Tilt relaxation: if true, relax the tilts and twist toward the ones
                                  !! the neighbours' means imply wherever the artificial viscosity acts.
  real :: dg_tilt_relax_frac      !< Tilt relaxation: its rate as a fraction of the artificial viscosity's
                                  !! jump decay rate on the cell's faces [nondim].
  real :: dg_tilt_relax_u_cut     !< Tilt relaxation: if positive, its rate is the speed law
                                  !! max(u_cut - u_credit*|u_n|, 0)/dx on each axis [L T-1 ~> m s-1].
  real :: dg_tilt_relax_u_credit  !< Tilt relaxation: the fraction of the transport's own removal
                                  !! speed that the speed law credits [nondim].
  logical :: dg_tilt_relax_twist  !< Tilt relaxation: if true, relax the twist as well as the tilts.
  logical :: dg_tilt_relax_free_edge !< Tilt relaxation: if true, also relax an axis on which the
                                  !! cell has no artificial-viscosity face on one side, so that its
                                  !! outer node has no pin and would otherwise drift.
  logical :: dg_tilt_relax_wide_nb !< Tilt relaxation: if true, a neighbour is usable when it holds
                                  !! ice, whatever its flotation state, and a cell the grounding
                                  !! line crosses is relaxed like any other.
  real, pointer, dimension(:,:) :: dg_tilt_relax_rate => NULL() !< Tilt relaxation: largest delivered
                                  !! rate over the three modes [T-1 ~> s-1].
  real, pointer, dimension(:,:) :: dg_tilt_relax_tend => NULL() !< Tilt relaxation: L2 cell norm of
                                  !! the thickness correction rate [Z T-1 ~> m s-1].
  real, pointer, dimension(:,:) :: dg_tilt_relax_alt_x => NULL() !< Tilt relaxation: alternating part
                                  !! of the x-tilt residual, where it is formed [Z ~> m].
  real, pointer, dimension(:,:) :: dg_tilt_relax_alt_y => NULL() !< Tilt relaxation: alternating part
                                  !! of the y-tilt residual, where it is formed [Z ~> m].
  real, pointer, dimension(:,:) :: dg_tilt_relax_alt_w => NULL() !< Tilt relaxation: alternating part
                                  !! of the twist residual, where it is formed [Z ~> m].
  real, pointer, dimension(:,:) :: dg_tilt_relax_state => NULL() !< Tilt relaxation: flotation state
                                  !! of each ice cell, 1 grounded, -1 floating, 0 crossed (skipped) [nondim].
  real :: dg_art_visc_r_hi        !< Artificial viscosity: surface jump over mean thickness at which the
                                  !! face coefficient reaches 1 [nondim].
  real :: dg_art_visc_kcell       !< Artificial viscosity: per-cell bound on dt times the summed face
                                  !! decay rates; SSP-RK2 needs < 2 [nondim].
  real :: dg_slow_idle_u_tiny     !< Artificial viscosity: stagnant-jump flag threshold on |u_face| [L T-1 ~> m s-1].
  real :: dg_slow_idle_eps_tiny   !< Artificial viscosity: stagnant-jump flag threshold on eps_e_face [T-1 ~> s-1].
  real :: dg_slow_idle_s_tol      !< Artificial viscosity: stagnant-jump flag threshold on |Delta h_eq| [Z ~> m].
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
  logical :: newton_divergence_rescue !< If true, monitor the nonlinear residual while Newton is
                                 !! active and, on divergence (residual NaN or exceeding
                                 !! newton_divergence_factor times its value at the Picard-to-
                                 !! Newton switch), restore the pre-Newton iterate, revert to
                                 !! Picard with a fresh outer-iteration budget, and reduce the
                                 !! switch threshold tenfold for the remainder of this solve.
  real :: newton_divergence_factor !< Factor on the nonlinear residual at the Picard-to-Newton
                                 !! switch above which the Newton iteration is declared
                                 !! divergent and rescued [nondim].
  integer :: newton_max_rescues  !< Maximum number of divergence rescues per velocity solve;
                                 !! once reached, Newton is disabled and the remainder of the
                                 !! solve runs pure Picard.
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
             id_dg_damp_tend = -1, id_dg_damp_gate = -1, &
             id_dg_damp_want = -1, id_dg_damp_got = -1, &
             id_dg_damp_ahat = -1, id_dg_damp_spread = -1, &
             id_dg_tilt_relax_rate = -1, id_dg_tilt_relax_tend = -1, &
             id_dg_tilt_relax_alt_x = -1, id_dg_tilt_relax_alt_y = -1, &
             id_dg_tilt_relax_alt_w = -1, id_dg_tilt_relax_state = -1, &
             id_ground_frac = -1, id_col_thick = -1, id_OD_av = -1, &
             id_f_ground_cell = -1, id_f_ground_node = -1, &
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
             id_h_source_rate = -1, &
             id_phi_x_FV = -1, id_phi_y_FV = -1, &
             id_dg_art_visc_coef_u = -1, id_dg_art_visc_coef_v = -1, &
             id_dg_art_visc_nu_u = -1, id_dg_art_visc_nu_v = -1, &
             id_dg_art_visc_cell_scale = -1, &
             id_dg_slow_idle_face_u = -1, id_dg_slow_idle_face_v = -1, &
             id_h_jump_face_u = -1, id_h_jump_face_v = -1, &
             id_s_jump_face_u = -1, id_s_jump_face_v = -1, &
             id_s_jump_face_u_rel = -1, id_s_jump_face_v_rel = -1, &
             id_h_jump_face_u_signed = -1, id_h_jump_face_v_signed = -1, &
             id_s_jump_face_u_signed = -1, id_s_jump_face_v_signed = -1, &
             id_un_face_u = -1, id_un_face_v = -1, &
             id_dg_eps_face_u = -1, id_dg_eps_face_v = -1
  real, pointer, dimension(:,:) :: dg_art_visc_coef_u => NULL() !< DG(1) art-visc coefficient on u-faces
                                                       !! after the per-cell cap [nondim].
  real, pointer, dimension(:,:) :: dg_art_visc_coef_v => NULL() !< DG(1) art-visc coefficient on v-faces
                                                       !! after the per-cell cap [nondim].
  real, pointer, dimension(:,:) :: dg_art_visc_nu_u => NULL() !< DG(1) art viscosity c_face*u_eff*dx_perp
                                                       !! on u-faces after the cap [L2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:) :: dg_art_visc_nu_v => NULL() !< DG(1) art viscosity c_face*u_eff*dx_perp
                                                       !! on v-faces after the cap [L2 T-1 ~> m2 s-1].
  real, pointer, dimension(:,:) :: dg_art_visc_cell_scale => NULL() !< Per-cell art-visc cap factor,
                                                       !! 1 where the budget is slack [nondim].
  real, pointer, dimension(:,:) :: dg_slow_idle_face_u => NULL() !< 1 on u-faces where |u_face|, eps_e_face
                                                       !! and |Delta h_eq| pass the stagnant-jump thresholds [nondim].
  real, pointer, dimension(:,:) :: dg_slow_idle_face_v => NULL() !< 1 on v-faces where |u_face|, eps_e_face
                                                       !! and |Delta h_eq| pass the stagnant-jump thresholds [nondim].
  real, pointer, dimension(:,:) :: phi_x_FV => NULL() !< Slope-limiter factor at each u-face
                                                       !! from ice_shelf_advect_thickness_x [nondim],
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
function slope_limiter(num, denom, limiter)
  real, intent(in)    :: num   !< The numerator of the ratio used in the slope limiter
  real, intent(in)    :: denom !< The denominator of the ratio used in the slope limiter
  integer, optional, intent(in) :: limiter !< The TVD limiter to use (LIMITER_VANLEER (default),
                                           !! LIMITER_SUPERBEE, LIMITER_MINMOD, or LIMITER_MC)
  real :: slope_limiter ! The slope limiter value, between 0 and 2 [nondim].
  real :: r  ! The ratio of num/denom [nondim]
  integer :: lim ! The selected limiter

  lim = LIMITER_VANLEER ; if (present(limiter)) lim = limiter

  if (denom == 0) then
    slope_limiter = 0
  elseif (num*denom <= 0) then  ! r <= 0: all TVD limiters return 0 at an extremum
    slope_limiter = 0
  else
    r = num/denom
    select case (lim)
      case (LIMITER_SUPERBEE)
        slope_limiter = max(min(2.0*r, 1.0), min(r, 2.0))
      case (LIMITER_MINMOD)
        slope_limiter = min(r, 1.0)
      case (LIMITER_MC)
        slope_limiter = min(min(2.0*r, 0.5*(1.0+r)), 2.0)
      case default  ! LIMITER_VANLEER
        slope_limiter = (r+abs(r))/(1+abs(r))
    end select
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
  logical :: DG_thickness   ! If true, DG(1) carries the thickness (USE_DG_THICKNESS)
  logical :: thickness_nodal ! If true, the thickness is read at nodes (INIT_ICE_THICKNESS_NODAL)
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
    allocate(CS%coef_prefactor_node(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%fB_node(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%area_node(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%OD_av(isd:ied,jsd:jed), source=0.0)
    allocate(CS%ground_frac(isd:ied,jsd:jed), source=0.0)
    ! Unconditional: the damper flags are not read yet.
    allocate(CS%dg_damp_ahat(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_damp_spread(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_damp_tend(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_damp_gate(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_damp_want(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_damp_got(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_tilt_relax_rate(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_tilt_relax_tend(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_tilt_relax_alt_x(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_tilt_relax_alt_y(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_tilt_relax_alt_w(isd:ied,jsd:jed), source=0.0)
    allocate(CS%dg_tilt_relax_state(isd:ied,jsd:jed), source=0.0)
    allocate(CS%f_ground_node(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%f_ground_cell(isd:ied,jsd:jed), source=0.0)
    allocate(CS%H_node(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%H_corner(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%fls_corner(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%corner_valid(IsdB:IedB,JsdB:JedB), source=.false.)
    allocate(CS%corner_wt(4,IsdB:IedB,JsdB:JedB), source=0.25)
    allocate(CS%taudx_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%taudy_shelf(IsdB:IedB,JsdB:JedB), source=0.0)
    allocate(CS%sx_shelf(isd:ied,jsd:jed), source=0.0)
    allocate(CS%sy_shelf(isd:ied,jsd:jed), source=0.0)
    allocate(CS%bed_elev(isd:ied,jsd:jed), source=0.0)
    allocate(CS%bed_node(IsdB:IedB,JsdB:JedB), source=0.0)
    ! Read early because they decide what is allocated and restarted. INIT_ICE_THICKNESS_NODAL
    ! also needs h_nodal, briefly, to derive the cell mean.
    call get_param(param_file, mdl, "USE_DG_THICKNESS", DG_thickness, &
                   default=.false., do_not_log=.true.)
    call get_param(param_file, mdl, "INIT_ICE_THICKNESS_NODAL", thickness_nodal, &
                   default=.false., do_not_log=.true.)
    if (DG_thickness .and. .not. thickness_nodal) call MOM_error(FATAL, &
        "MOM_ice_shelf_dynamics: USE_DG_THICKNESS initializes the DG(1) thickness from nodal "//&
        "values; set INIT_ICE_THICKNESS_NODAL=True.")
    if (DG_thickness .or. thickness_nodal) &
      allocate(CS%h_nodal(isd:ied,jsd:jed,1:2,1:2), source=0.0)
    if (DG_thickness) then
      allocate(CS%Minv_nodal(isd:ied,jsd:jed,1:2,1:2,1:2,1:2), source=0.0)
      allocate(CS%cell_mean_w(isd:ied,jsd:jed,1:2,1:2), source=0.0)
    endif
    allocate(CS%h_source_rate(isd:ied,jsd:jed), source=0.0)
    allocate(CS%h_source_rate_bmb(isd:ied,jsd:jed), source=0.0)
    allocate(CS%xi_basal(isd:ied,jsd:jed,1:2,1:2), source=1.0)
    allocate(CS%h_source_rate_last(isd:ied,jsd:jed), source=0.0)
    allocate(CS%phi_x_FV(IsdB:IedB,jsd:jed), source=1.0)
    allocate(CS%phi_y_FV(isd:ied,JsdB:JedB), source=1.0)
    allocate(CS%dg_art_visc_coef_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_art_visc_coef_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%dg_art_visc_nu_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_art_visc_nu_v(isd:ied,JsdB:JedB), source=0.0)
    allocate(CS%dg_art_visc_cell_scale(isd:ied,jsd:jed), source=1.0)
    allocate(CS%dg_slow_idle_face_u(IsdB:IedB,jsd:jed), source=0.0)
    allocate(CS%dg_slow_idle_face_v(isd:ied,JsdB:JedB), source=0.0)
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
    ! Under DG, h_shelf is derived from h_nodal, so only the corners need restarting.
    if (DG_thickness) then
      call register_restart_field(CS%h_nodal(:,:,1,1), "h_nodal_SW_DG", .true., restart_CS, &
                                  "DG(1) nodal Q1 thickness, SW corner", "m", conversion=US%Z_to_m)
      call register_restart_field(CS%h_nodal(:,:,2,1), "h_nodal_SE_DG", .true., restart_CS, &
                                  "DG(1) nodal Q1 thickness, SE corner", "m", conversion=US%Z_to_m)
      call register_restart_field(CS%h_nodal(:,:,1,2), "h_nodal_NW_DG", .true., restart_CS, &
                                  "DG(1) nodal Q1 thickness, NW corner", "m", conversion=US%Z_to_m)
      call register_restart_field(CS%h_nodal(:,:,2,2), "h_nodal_NE_DG", .true., restart_CS, &
                                  "DG(1) nodal Q1 thickness, NE corner", "m", conversion=US%Z_to_m)
    endif
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
  character(len=200) :: IS_energyfile  ! The name of the energy file.
  character(len=32) :: filename_appendix = '' ! FMS appendix to filename for ensemble runs
  character(len=16) :: inner_solver_str ! The type of inner solver to use for the SSA
  character(len=16) :: gl_subgrid_scheme_str ! Grounding-line subgrid quadrature scheme string
  character(len=16) :: adv_limiter_str ! Thickness-advection TVD slope-limiter choice string
  character(len=16) :: melt_glp_str    ! Ice-only prescribed basal melt grounding-line scheme string
  character(len=16) :: cism_friction_str ! CISM-style basal friction mode string
  character(len=16) :: cism_taud_str     ! CISM-style driving stress mode string
  logical :: solo_ice_sheet   ! True if this is an ice-only (solo ice sheet) run

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
    call get_param(param_file, mdl, "GROUNDING_LINE_SUBGRID_SCHEME", gl_subgrid_scheme_str, &
                 "Quadrature scheme for the sub-cell grounding-line integration when "//&
                 "GROUNDING_LINE_INTERPOLATE is true. 'SEP3' samples each cell with a uniform "//&
                 "NxN sub-grid (GROUNDING_LINE_INTERP_SUBGRID_N) and a per-point flotation "//&
                 "test. 'SEP2' splits each grounding-line cell geometrically into grounded and "//&
                 "floating sub-elements and integrates each side exactly on its own quadrature "//&
                 "(Seroussi et al. 2014, extended to quadrilateral elements); it applies to the "//&
                 "basal friction, the DG driving stress, and the grounded fraction, and ignores "//&
                 "GROUNDING_LINE_INTERP_SUBGRID_N.", &
                 default="SEP3", do_not_log=.not.CS%GL_regularize)
    select case (trim(gl_subgrid_scheme_str))
      case ("SEP3") ; CS%use_sep2 = .false.
      case ("SEP2") ; CS%use_sep2 = .true.
      case default  ; call MOM_error(FATAL, "MOM_ice_shelf_dynamics: "//&
                        "GROUNDING_LINE_SUBGRID_SCHEME must be 'SEP3' or 'SEP2'.")
    end select
    call get_param(param_file, mdl, "GROUNDING_LINE_COUPLE", CS%GL_couple, &
                 "If true, let the floatation condition be determined by "//&
                 "ocean column thickness. This means that update_OD_ffrac "//&
                 "will be called.  GL_REGULARIZE and GL_COUPLE are exclusive.", &
                 default=.false., do_not_log=CS%GL_regularize)
    if (CS%GL_regularize) CS%GL_couple = .false.
    if (present(solo_ice_sheet_in)) then
      if (solo_ice_sheet_in) CS%GL_couple = .false.
    endif
    if (CS%GL_regularize .and. (CS%n_sub_regularize == 0)) then
      if (CS%use_sep2) then
        ! SEP2 does not use the uniform sub-grid; keep Phisub minimally allocated.
        CS%n_sub_regularize = 1
      else
        call MOM_error (FATAL, &
          "GROUNDING_LINE_INTERP_SUBGRID_N must be a positive integer if GL regularization is used")
      endif
    endif
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
    call get_param(param_file, mdl, "FV_GL_ONE_SIDED_TAUD", CS%FV_GL_one_sided, &
                 "If true, the finite-volume (non-DG) driving stress is evaluated with "//&
                 "one-sided differences in the grounded and floating cells that border the "//&
                 "grounding line, following Cornford et al. (2013) eqs 27-29, instead of a "//&
                 "centered difference that straddles the grounding line.", &
                 default=.false.)
    call get_param(param_file, mdl, "CISM_FRICTION", cism_friction_str, &
                 "How the basal friction treats the grounding line, following Leguy, Lipscomb "//&
                 "& Asay-Davis (2021, The Cryosphere 15:3229-3253) and CISM. 'off' leaves the "//&
                 "friction to the geometric sub-cell integration selected by "//&
                 "GROUNDING_LINE_SUBGRID_SCHEME. 'local' and 'integrate' both scale the friction "//&
                 "by the analytic nodal grounded fraction of the quadrant parameterization "//&
                 "(sec. 2.2), which is the bilinear-flotation area integral and so is rotation-"//&
                 "consistent and needs no sub-cell sampling. They differ in the assembly: "//&
                 "'local' is the nodal diagonal of CISM HO_ASSEMBLE_BETA_LOCAL, where the drag "//&
                 "at each node is beta(node)*area_node*u(node) with no element integration or "//&
                 "neighbor coupling, and 'integrate' is the consistent element-quadrature mass. "//&
                 "'local' is the CISM-faithful choice and co-locates with CISM_TAUD='local'.", &
                 default="off")
    select case (trim(cism_friction_str))
      case ("off")       ; CS%cism_friction = CISM_OFF
      case ("local")     ; CS%cism_friction = CISM_LOCAL
      case ("integrate") ; CS%cism_friction = CISM_INTEGRATE
      case default ; call MOM_error(FATAL, "MOM_ice_shelf_dynamics: CISM_FRICTION must be "//&
                 "'off', 'local' or 'integrate', but got '"//trim(cism_friction_str)//"'.")
    end select
    CS%gl_quad_friction     = (CS%cism_friction /= CISM_OFF)
    CS%local_basal_friction = (CS%cism_friction == CISM_LOCAL)
    call get_param(param_file, mdl, "USE_DG_THICKNESS", CS%use_DG_thickness, &
                 "If true, represent ice thickness as discontinuous bilinear (DG(1)) corner "//&
                 "values per cell, advected with SSP-RK2.", &
                 default=.false.)
    call get_param(param_file, mdl, "INIT_ICE_BED_NODAL", CS%init_bed_nodal, &
                 "If true, read BED_TOPO_VARNAME at B-grid nodes and set the cell bed to the area-weighted "//&
                 "mean of their bilinear interpolant. Required by USE_DG_THICKNESS; optional otherwise.", &
                 default=.false.)
    if (CS%use_DG_thickness .and. .not. CS%init_bed_nodal) call MOM_error(FATAL, &
        "MOM_ice_shelf_dynamics: USE_DG_THICKNESS needs the bed at B-grid nodes to evaluate "//&
        "grad(b) inside an element; set INIT_ICE_BED_NODAL=True.")

    ! Prescribed basal melt for the ice-only driver.
    solo_ice_sheet = .false.
    if (present(solo_ice_sheet_in)) solo_ice_sheet = solo_ice_sheet_in
    CS%dg_basal_source_sem2 = .false.
    call get_param(param_file, mdl, "ICE_ONLY_BASAL_MELT", CS%ice_only_basal_melt, &
                 "If true, the ice-only driver prescribes a basal melt (Seroussi & Morlighem 2018, eq 4; "//&
                 "Leguy et al. 2021 eq 18): 0 at an ice-base depth of 50 m rising linearly to 30 m yr-1 at >=500 m.", &
                 default=.false., do_not_log=.not.solo_ice_sheet)
    call get_param(param_file, mdl, "ICE_ONLY_BASAL_MELT_GLP", melt_glp_str, &
                 "How the ice-only basal melt is applied in grounding-line cells (Leguy et al. "//&
                 "2021; Seroussi & Morlighem 2018). 'FMP': full rate. 'FCMP': full rate where the "//&
                 "cell centre floats. 'PMP': rate times the floating fraction. 'NMP': none in "//&
                 "partly grounded cells. 'SEM2': PMP's cell total, distributed by nodal floating "//&
                 "fraction; requires USE_DG_THICKNESS and GROUNDING_LINE_INTERPOLATE or CISM_FRICTION.", &
                 default="FMP", do_not_log=.not.CS%ice_only_basal_melt)
    select case (trim(melt_glp_str))
      case ("FMP")  ; CS%ice_only_melt_glp = MELT_GLP_FMP
      case ("FCMP") ; CS%ice_only_melt_glp = MELT_GLP_FCMP
      case ("PMP")  ; CS%ice_only_melt_glp = MELT_GLP_PMP
      case ("NMP")  ; CS%ice_only_melt_glp = MELT_GLP_NMP
      case ("SEM2") ; CS%ice_only_melt_glp = MELT_GLP_SEM2
      case default  ; call MOM_error(FATAL, "MOM_ice_shelf_dynamics: "//&
                        "ICE_ONLY_BASAL_MELT_GLP must be one of 'FMP', 'FCMP', 'PMP', 'NMP' or "//&
                        "'SEM2', but got '"//trim(melt_glp_str)//"'.")
    end select
    CS%dg_basal_source_sem2 = (CS%ice_only_melt_glp == MELT_GLP_SEM2)
    call get_param(param_file, mdl, "ICE_ONLY_BASAL_MELT_SCALE", CS%ice_only_melt_scale, &
                 "Factor on the prescribed ice-only basal melt rate; 5 gives the high-melt experiments of "//&
                 "Leguy et al. (2021).", &
                 units="nondim", default=1.0, do_not_log=.not.CS%ice_only_basal_melt)
    if (CS%ice_only_melt_scale < 0.0) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: ICE_ONLY_BASAL_MELT_SCALE must be non-negative.")
    if (CS%ice_only_basal_melt .and. .not.solo_ice_sheet) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: ICE_ONLY_BASAL_MELT is only meaningful for the ice-only "//&
                 "driver; in a coupled run the basal melt rate is supplied by the ocean.")
    if (CS%ice_only_melt_glp == MELT_GLP_SEM2) then
      if (.not.CS%use_DG_thickness) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: ICE_ONLY_BASAL_MELT_GLP = 'SEM2' requires USE_DG_THICKNESS.")
      if (.not.(CS%GL_regularize .or. CS%gl_quad_friction)) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: ICE_ONLY_BASAL_MELT_GLP = 'SEM2' requires "//&
                 "GROUNDING_LINE_INTERPOLATE or CISM_FRICTION.")
    endif
    call get_param(param_file, mdl, "CISM_TAUD", cism_taud_str, &
                 "How the finite-volume (non-DG) driving stress treats the grounding line, "//&
                 "following Lipscomb et al. (2019, Geosci. Model Dev. 12:387-424) and CISM. "//&
                 "'off' uses the cell-centroid centered difference of the surface. 'local' and "//&
                 "'integrate' both evaluate the surface gradient directly at B-grid nodes from "//&
                 "the four surrounding cell centers (eq. 14, with their 'option 3' ice-margin "//&
                 "treatment), a compact stencil that is less smeared across the grounding line. "//&
                 "They differ in the assembly: 'local' is the mass-lumped method of A4, where "//&
                 "each node's driving stress uses the surface slope at that node alone over its "//&
                 "nodal control mass (CISM HO_ASSEMBLE_TAUD_LOCAL), and 'integrate' is the "//&
                 "consistent element-quadrature assembly. 'local' is the CISM-faithful choice "//&
                 "and co-locates with CISM_FRICTION='local'. Honors MAX_SURFACE_SLOPE; not "//&
                 "'off' is mutually exclusive with FV_GL_ONE_SIDED_TAUD.", &
                 default="off")
    select case (trim(cism_taud_str))
      case ("off")       ; CS%cism_taud = CISM_OFF
      case ("local")     ; CS%cism_taud = CISM_LOCAL
      case ("integrate") ; CS%cism_taud = CISM_INTEGRATE
      case default ; call MOM_error(FATAL, "MOM_ice_shelf_dynamics: CISM_TAUD must be "//&
                 "'off', 'local' or 'integrate', but got '"//trim(cism_taud_str)//"'.")
    end select
    CS%fv_taud_vertex_grad  = (CS%cism_taud /= CISM_OFF)
    CS%local_fv_taud_vertex = (CS%cism_taud == CISM_LOCAL)
    if (CS%fv_taud_vertex_grad .and. CS%FV_GL_one_sided) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: CISM_TAUD replaces the cell-centroid surface slope with "//&
                 "a nodal gradient, which has no one-sided analog; it cannot be used with "//&
                 "FV_GL_ONE_SIDED_TAUD.")

    ! Sub-element grounding line for the FV (non-DG) path: one flotation field (CS%fls_corner) on one
    ! partition drives the grounding-line location, the basal friction, the Coulomb effective pressure,
    ! and the driving stress. Both parameters require GROUNDING_LINE_INTERPOLATE, which is what
    ! allocates Phisub and enables compute_ground_frac and the sub-element quadrature dispatch; without
    ! it there is no sub-cell partition for either term to integrate over.
    call get_param(param_file, mdl, "I2N_FRICTION", CS%i2n_friction, &
                 "If true, integrate the finite-volume (non-DG) basal friction and the Coulomb "//&
                 "effective pressure over the sub-element grounding-line partition selected by "//&
                 "GROUNDING_LINE_SUBGRID_SCHEME, using the corner thickness and flotation deficit "//&
                 "obtained by interpolating the cell-centered fields with dual-cell Lagrange weights "//&
                 "over ice-covered cells (plus ice-free land lying below the ice). This replaces the "//&
                 "corner-thickness-with-cell-constant-bed flotation field, so the grounding line seen "//&
                 "by the friction is the same one seen by I2N_TAUD. The effective pressure "//&
                 "is evaluated at each grounded quadrature point as rho_ocean*g*min(fls, r*H). "//&
                 "Requires GROUNDING_LINE_INTERPOLATE=True.", &
                 default=.false.)
    call get_param(param_file, mdl, "I2N_TAUD", CS%i2n_taud, &
                 "If true, integrate the finite-volume (non-DG) driving stress over the same "//&
                 "sub-element grounding-line partition used by I2N_FRICTION. The surface "//&
                 "elevation is reconstructed as S = (1-r)*H + max(fls,0) from the same two corner "//&
                 "fields, so the slope kink lies exactly on the partition's grounding line and every "//&
                 "quadrature point takes one side of it; the cell-mean thickness multiplies the "//&
                 "resulting slope. Requires GROUNDING_LINE_INTERPOLATE=True.", &
                 default=.false.)
    if ((CS%i2n_friction .or. CS%i2n_taud) .and. .not. CS%GL_regularize) &
      call MOM_error(FATAL, "MOM_ice_shelf_dynamics: I2N_FRICTION and I2N_TAUD "//&
                 "integrate over the sub-cell grounding-line partition and require "//&
                 "GROUNDING_LINE_INTERPOLATE=True.")
    ! The interpolate-to-nodes and CISM-style schemes are two different sources of the
    ! grounded fraction, and 'local' is additionally a nodal diagonal on the dual cell
    ! where I2N assembles over the primal element with Q1 weighting, so mixing them
    ! mismatches the control volumes at the grounding line.
    if (CS%i2n_friction .and. (CS%cism_friction /= CISM_OFF)) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: I2N_FRICTION and CISM_FRICTION are two different "//&
                 "grounding-line treatments of the basal friction and are mutually exclusive; "//&
                 "set CISM_FRICTION='off'.")
    if (CS%i2n_taud .and. (CS%cism_taud /= CISM_OFF)) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: I2N_TAUD and CISM_TAUD are two different "//&
                 "grounding-line treatments of the driving stress and are mutually exclusive; "//&
                 "set CISM_TAUD='off'.")
    if (CS%i2n_taud .and. .not. CS%use_sep2) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: I2N_TAUD integrates the surface-slope kink on "//&
                 "the geometric sub-element partition and currently requires "//&
                 "GROUNDING_LINE_SUBGRID_SCHEME='SEP2'. I2N_FRICTION supports both "//&
                 "schemes.")
    if (CS%i2n_taud .and. CS%FV_GL_one_sided) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: I2N_TAUD replaces the cell-centroid surface slope "//&
                 "with a sub-element reconstruction, which has no one-sided analog; it cannot be used "//&
                 "with FV_GL_ONE_SIDED_TAUD.")

    ! CISM_* only read cell means, so they work under DG. I2N would discard the nodal thickness.
    if (CS%use_DG_thickness .and. CS%i2n_friction) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: I2N_FRICTION reconstructs a corner thickness from the "//&
                 "cell means, which would discard the DG(1) nodal thickness; it cannot be used "//&
                 "with USE_DG_THICKNESS.")
    if (CS%use_DG_thickness .and. CS%i2n_taud) call MOM_error(FATAL, &
                 "MOM_ice_shelf_dynamics: I2N_TAUD reconstructs a corner surface from the cell "//&
                 "means, which would discard the DG(1) nodal thickness; it cannot be used with "//&
                 "USE_DG_THICKNESS.")
    if (CS%use_DG_thickness .and. (CS%cism_friction /= CISM_OFF)) call MOM_error(WARNING, &
                 "MOM_ice_shelf_dynamics: USE_DG_THICKNESS with CISM_FRICTION takes the grounded "//&
                 "fraction from the cell-mean quadrant parameterization rather than the DG(1) "//&
                 "nodal thickness; DG normally supplies the basal friction on its own basis.")
    if (CS%use_DG_thickness .and. (CS%cism_taud /= CISM_OFF)) call MOM_error(WARNING, &
                 "MOM_ice_shelf_dynamics: USE_DG_THICKNESS with CISM_TAUD routes the driving "//&
                 "stress through the finite-volume path on the cell-mean thickness, discarding "//&
                 "the DG(1) sub-cell surface slope; DG normally supplies the driving stress on "//&
                 "its own basis.")
    call get_param(param_file, mdl, "ICE_SHELF_ADVECT_LIMITER", adv_limiter_str, &
                 "The TVD slope limiter used for the finite-volume ice thickness advection in "//&
                 "ice_shelf_advect_thickness_x/y. VAN_LEER is the original scheme; SUPERBEE is "//&
                 "the least diffusive (matches the STREAMICE default) but compressive; MINMOD is "//&
                 "the most diffusive; MC is intermediate.", default="VAN_LEER")
    select case (trim(adv_limiter_str))
      case ("VAN_LEER") ; CS%adv_thickness_limiter = LIMITER_VANLEER
      case ("SUPERBEE") ; CS%adv_thickness_limiter = LIMITER_SUPERBEE
      case ("MINMOD")   ; CS%adv_thickness_limiter = LIMITER_MINMOD
      case ("MC")       ; CS%adv_thickness_limiter = LIMITER_MC
      case default ; call MOM_error(FATAL, "ICE_SHELF_ADVECT_LIMITER = "//trim(adv_limiter_str)//&
                 " is invalid; use VAN_LEER, SUPERBEE, MINMOD, or MC.")
    end select
    call get_param(param_file, mdl, "ICE_SHELF_ADVECT_CFL_WEIGHT", CS%adv_cfl_weight, &
                 "If true, weight the ice thickness-advection slope reconstruction by the "//&
                 "Lax-Wendroff (1-CFL) factor, as in STREAMICE and standard flux-form TVD "//&
                 "schemes, giving a time-accurate 2nd-order flux. If false, the slope uses the "//&
                 "full spatial reconstruction (the original behavior).", default=.false.)
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
    call get_param(param_file, mdl, "NEWTON_DIVERGENCE_RESCUE", CS%newton_divergence_rescue, &
                "If true, monitor the nonlinear residual while Newton iterations are active "//&
                "and, if it becomes NaN or exceeds NEWTON_DIVERGENCE_FACTOR times its value "//&
                "at the Picard-to-Newton switch, restore the pre-Newton velocity iterate, "//&
                "revert to Picard iterations with a fresh outer-iteration budget, and reduce "//&
                "the Picard-to-Newton switch threshold by a factor of 10 for the remainder "//&
                "of this velocity solve (the configured NEWTON_AFTER_TOLERANCE is restored "//&
                "at the next solve). At most NEWTON_DIVERGENCE_MAX_RESCUES rescues are "//&
                "attempted per solve, after which Newton is disabled and the solve "//&
                "completes as pure Picard. No effect when NEWTON_AFTER_TOLERANCE <= 0.", &
                default=.false.)
    call get_param(param_file, mdl, "NEWTON_DIVERGENCE_FACTOR", CS%newton_divergence_factor, &
                "Factor on the nonlinear residual at the Picard-to-Newton switch above "//&
                "which the Newton iteration is declared divergent and rescued.", &
                units="nondim", default=10.0, do_not_log=.not.CS%newton_divergence_rescue)
    call get_param(param_file, mdl, "NEWTON_DIVERGENCE_MAX_RESCUES", CS%newton_max_rescues, &
                "Maximum number of Newton divergence rescues per velocity solve. Once "//&
                "reached, Newton is disabled (the working switch threshold is set to "//&
                "zero) and the remainder of the solve runs pure Picard.", &
                default=2, do_not_log=.not.CS%newton_divergence_rescue)
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
    call read_DG_params(param_file, mdl, CS, US)
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

    ! Per-cell nodal DG(1) metric tables (Minv_nodal, cell_mean_w).
    if (CS%use_DG_thickness) call init_nodal_DG_metric(CS, G)

    if (CS%GL_regularize) then
      allocate(CS%Phisub(2,2,CS%n_sub_regularize,CS%n_sub_regularize,2,2), source=0.0)
      call bilinear_shape_functions_subgrid(CS%Phisub, CS%n_sub_regularize)
    endif

    ! Dual-cell Lagrange weights for the FV sub-element corner fields; grid-only, so once at init.
    if (CS%i2n_friction .or. CS%i2n_taud) &
      call build_corner_lagrange_weights(CS, G)

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
      call initialize_bed_node_from_file(CS%bed_node, CS%bed_elev, G, US, param_file)
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
                  G, US, param_file, skip_bed=CS%init_bed_nodal)
      call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE, complete=.true.)
      call pass_var(CS%ground_frac, G%domain, complete=.true.)
      if (CS%init_bed_nodal) then
        ! Reads bed_node from file and derives bed_elev (both halo-updated inside).
        call initialize_bed_node_from_file(CS%bed_node, CS%bed_elev, G, US, param_file)
      else
        call pass_var(CS%bed_elev, G%domain, complete=.true.)
      endif
      if (CS%use_DG_thickness) then
        ! DG(1) cold start from the nodal thickness read by initialize_ice_thickness.
        call pass_var(ISS%h_shelf, G%domain)
        call nodal_positivity_limit(CS, G, ISS)
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
    if ((CS%dg_tilt_damp .or. CS%dg_twist_damp) .and. associated(CS%dg_damp_tend)) then
      CS%id_dg_damp_tend = register_diag_field('ice_shelf_model','dg_mode_damp_tend', &
         CS%diag%axesT1, Time, 'L2 cell norm of the DG(1) tilt-mode-damper thickness correction rate', &
         'm s-1', conversion=US%Z_to_m*US%s_to_T)
      CS%id_dg_damp_gate = register_diag_field('ice_shelf_model','dg_mode_damp_gate', &
         CS%diag%axesT1, Time, 'largest DG(1) tilt-mode-damper gate over its three modes, 0 to 1', 'none')
      CS%id_dg_damp_want = register_diag_field('ice_shelf_model','dg_mode_damp_rate_want', &
         CS%diag%axesT1, Time, 'DG(1) tilt-mode-damper rate requested, summed over modes, before '//&
         'the transport credit', 's-1', conversion=US%s_to_T)
      CS%id_dg_damp_got = register_diag_field('ice_shelf_model','dg_mode_damp_rate_got', &
         CS%diag%axesT1, Time, 'DG(1) tilt-mode-damper rate delivered, summed over modes', &
         's-1', conversion=US%s_to_T)
      CS%id_dg_damp_ahat = register_diag_field('ice_shelf_model','dg_mode_damp_ahat', &
         CS%diag%axesT1, Time, 'L2 cell norm of the DG(1) mode-damper detector over its three '//&
         'modes', 'm', conversion=US%Z_to_m)
      CS%id_dg_damp_spread = register_diag_field('ice_shelf_model','dg_mode_damp_spread', &
         CS%diag%axesT1, Time, 'largest spread of the DG(1) mode-damper detector readings over '//&
         'its three modes; zero unless the agreement detector is selected', 'm', &
         conversion=US%Z_to_m)
    endif
    if (CS%dg_tilt_relax) then
      CS%id_dg_tilt_relax_rate = register_diag_field('ice_shelf_model','dg_tilt_relax_rate', &
         CS%diag%axesT1, Time, 'largest DG(1) tilt-relaxation rate delivered over the three modes', &
         's-1', conversion=US%s_to_T)
      CS%id_dg_tilt_relax_tend = register_diag_field('ice_shelf_model','dg_tilt_relax_tend', &
         CS%diag%axesT1, Time, 'L2 cell norm of the DG(1) tilt-relaxation thickness correction rate', &
         'm s-1', conversion=US%Z_to_m*US%s_to_T)
      CS%id_dg_tilt_relax_alt_x = register_diag_field('ice_shelf_model','dg_tilt_relax_alt_x', &
         CS%diag%axesT1, Time, 'alternating part of the DG(1) x-tilt residual seen by the tilt '//&
         'relaxation, in surface elevation for an unpinned cell the grounding line crosses; '//&
         'zero where it is not formed', 'm', conversion=US%Z_to_m)
      CS%id_dg_tilt_relax_alt_y = register_diag_field('ice_shelf_model','dg_tilt_relax_alt_y', &
         CS%diag%axesT1, Time, 'alternating part of the DG(1) y-tilt residual seen by the tilt '//&
         'relaxation, in surface elevation for an unpinned cell the grounding line crosses; '//&
         'zero where it is not formed', 'm', conversion=US%Z_to_m)
      CS%id_dg_tilt_relax_alt_w = register_diag_field('ice_shelf_model','dg_tilt_relax_alt_w', &
         CS%diag%axesT1, Time, 'alternating part of the DG(1) twist residual seen by the tilt '//&
         'relaxation; zero where it is not formed', 'm', conversion=US%Z_to_m)
      CS%id_dg_tilt_relax_state = register_diag_field('ice_shelf_model','dg_tilt_relax_state', &
         CS%diag%axesT1, Time, 'flotation state used by the DG(1) tilt relaxation: 1 grounded, '//&
         '-1 floating, 0 partly grounded (skipped, or relaxed in surface form along an '//&
         'unpinned axis) or no ice', 'nondim')
    endif

    CS%id_ground_frac = register_diag_field('ice_shelf_model','ice_ground_frac',CS%diag%axesT1, Time, &
       'fraction of cell that is grounded; under GL_regularize this is the fraction of '//&
       'sub-cell quadrature points whose draft sits below the bed', 'none')
    CS%id_f_ground_cell = register_diag_field('ice_shelf_model','f_ground_cell',CS%diag%axesT1, Time, &
      'analytic grounded ice fraction at cell centers from the quadrant grounding-line '//&
      'parameterization (Leguy et al. 2021); nonzero only when CISM_FRICTION is set', 'none')
    CS%id_f_ground_node = register_diag_field('ice_shelf_model','f_ground_node',CS%diag%axesB1, Time, &
      'analytic grounded ice fraction at B-grid nodes from the quadrant grounding-line '//&
      'parameterization (Leguy et al. 2021); multiplies basal friction under '//&
      'CISM_FRICTION', 'none')
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
         'Bed elevation at B-grid nodes', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_SW = register_diag_field('ice_shelf_model','h_nodal_SW',CS%diag%axesT1, Time, &
         'DG(1) thickness at SW cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_SE = register_diag_field('ice_shelf_model','h_nodal_SE',CS%diag%axesT1, Time, &
         'DG(1) thickness at SE cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_NW = register_diag_field('ice_shelf_model','h_nodal_NW',CS%diag%axesT1, Time, &
         'DG(1) thickness at NW cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_nodal_NE = register_diag_field('ice_shelf_model','h_nodal_NE',CS%diag%axesT1, Time, &
         'DG(1) thickness at NE cell corner', 'm', conversion=US%Z_to_m)
      CS%id_h_jump_node = register_diag_field('ice_shelf_model','h_jump_node',CS%diag%axesB1, Time, &
         'DG(1) max minus min corner thickness of the ice cells at a B-grid node', &
         'm', conversion=US%Z_to_m)
      CS%id_h_jump_node_rel = register_diag_field('ice_shelf_model','h_jump_node_rel',CS%diag%axesB1, Time, &
         'h_jump_node over the mean thickness of those cells', 'nondim')
      CS%id_h_node_max = register_diag_field('ice_shelf_model','h_node_max',CS%diag%axesB1, Time, &
         'DG(1) max corner thickness of the ice cells at a B-grid node', 'm', conversion=US%Z_to_m)
      CS%id_h_node_min = register_diag_field('ice_shelf_model','h_node_min',CS%diag%axesB1, Time, &
         'DG(1) min corner thickness of the ice cells at a B-grid node', 'm', conversion=US%Z_to_m)
      CS%id_h_source_rate = register_diag_field('ice_shelf_model','h_source_rate',CS%diag%axesT1, Time, &
         'Cell-mean thickness source (basal + surface) used by the last DG advection step', &
         'm s-1', conversion=US%Z_to_m*US%s_to_T)
      CS%id_dg_art_visc_coef_u = register_diag_field('ice_shelf_model','dg_art_visc_coef_u', &
         CS%diag%axesCu1, Time, &
         'DG(1) artificial-viscosity coefficient on u-faces, after the per-cell cap.'//&
         'Ranges from 0 to DG1_ART_VISC_C_MAX]', 'nondim')
      CS%id_dg_art_visc_coef_v = register_diag_field('ice_shelf_model','dg_art_visc_coef_v', &
         CS%diag%axesCv1, Time, &
         'DG(1) artificial-viscosity coefficient on v-faces, after the per-cell cap' //&
         'Ranges from 0 to DG1_ART_VISC_C_MAX]', 'nondim')
      CS%id_dg_art_visc_nu_u = register_diag_field('ice_shelf_model','dg_art_visc_nu_u', &
         CS%diag%axesCu1, Time, &
         'DG(1) artificial viscosity c_face*u_eff*dx_perp on u-faces, after the per-cell cap', &
         'm2 s-1', conversion=US%L_T_to_m_s*US%L_to_m)
      CS%id_dg_art_visc_nu_v = register_diag_field('ice_shelf_model','dg_art_visc_nu_v', &
         CS%diag%axesCv1, Time, &
         'DG(1) artificial viscosity c_face*u_eff*dx_perp on v-faces, after the per-cell cap', &
         'm2 s-1', conversion=US%L_T_to_m_s*US%L_to_m)
      CS%id_dg_art_visc_cell_scale = register_diag_field('ice_shelf_model', &
         'dg_art_visc_cell_scale', CS%diag%axesT1, Time, &
         'DG(1) artificial-viscosity cap factor: 1 where the DG1_ART_VISC_KCELL cap is inactive, '//&
         '< 1 cap engaged and face coefficients scaled by this factor', 'nondim')
      CS%id_dg_slow_idle_face_u = register_diag_field('ice_shelf_model','dg_slow_idle_face_u', &
         CS%diag%axesCu1, Time, &
         '1 on u-faces with a jump where speed and strain rate — and their damping — are negligible','nondim')
      CS%id_dg_slow_idle_face_v = register_diag_field('ice_shelf_model','dg_slow_idle_face_v', &
         CS%diag%axesCv1, Time, &
         '1 on v-faces with a jump where speed and strain rate — and their damping — are negligible','nondim')
      CS%id_h_jump_face_u = register_diag_field('ice_shelf_model','h_jump_face_u', &
         CS%diag%axesCu1, Time, &
         'DG(1) thickness jump across u-faces, max |[h]| over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_h_jump_face_v = register_diag_field('ice_shelf_model','h_jump_face_v', &
         CS%diag%axesCv1, Time, &
         'DG(1) thickness jump across v-faces, max |[h]| over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_u = register_diag_field('ice_shelf_model','s_jump_face_u', &
         CS%diag%axesCu1, Time, &
         'DG(1) surface jump across u-faces, max |[s]| over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_v = register_diag_field('ice_shelf_model','s_jump_face_v', &
         CS%diag%axesCv1, Time, &
         'DG(1) surface jump across v-faces, max |[s]| over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_u_rel = register_diag_field('ice_shelf_model','s_jump_face_u_rel', &
         CS%diag%axesCu1, Time, &
         's_jump_face_u (max |[s]| over the 2 face nodes) over the mean thickness of the two cells', &
         'nondim')
      CS%id_s_jump_face_v_rel = register_diag_field('ice_shelf_model','s_jump_face_v_rel', &
         CS%diag%axesCv1, Time, &
         's_jump_face_v (max |[s]| over the 2 face nodes) over the mean thickness of the two cells', 'nondim')
      CS%id_h_jump_face_u_signed = register_diag_field('ice_shelf_model','h_jump_face_u_signed', &
         CS%diag%axesCu1, Time, &
         'DG(1) thickness jump east minus west on u-faces, mean over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_h_jump_face_v_signed = register_diag_field('ice_shelf_model','h_jump_face_v_signed', &
         CS%diag%axesCv1, Time, &
         'DG(1) thickness jump north minus south on v-faces, mean over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_u_signed = register_diag_field('ice_shelf_model','s_jump_face_u_signed', &
         CS%diag%axesCu1, Time, &
         'DG(1) surface jump east minus west on u-faces, mean over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_s_jump_face_v_signed = register_diag_field('ice_shelf_model','s_jump_face_v_signed', &
         CS%diag%axesCv1, Time, &
         'DG(1) surface jump north minus south on v-faces, mean over the 2 face nodes', &
         'm', conversion=US%Z_to_m)
      CS%id_un_face_u = register_diag_field('ice_shelf_model','un_face_u', &
         CS%diag%axesCu1, Time, &
         '|u| on u-faces, mean of the 2 face nodes', 'm s-1', conversion=US%L_T_to_m_s)
      CS%id_un_face_v = register_diag_field('ice_shelf_model','un_face_v', &
         CS%diag%axesCv1, Time, &
         '|v| on v-faces, mean of the 2 face nodes', 'm s-1', conversion=US%L_T_to_m_s)
      CS%id_dg_eps_face_u = register_diag_field('ice_shelf_model','dg_eps_face_u', &
         CS%diag%axesCu1, Time, &
         'SSA effective strain rate at u-faces, as used by the DG(1) artificial viscosity', &
         'yr-1', conversion=365.0*86400.0*US%s_to_T)
      CS%id_dg_eps_face_v = register_diag_field('ice_shelf_model','dg_eps_face_v', &
         CS%diag%axesCv1, Time, &
         'SSA effective strain rate at v-faces, as used by the DG(1) artificial viscosity', &
         'yr-1', conversion=365.0*86400.0*US%s_to_T)
    else
      CS%id_phi_x_FV = register_diag_field('ice_shelf_model','phi_x_FV',CS%diag%axesCu1, Time, &
         'Slope-limiter factor at each u-face from ice_shelf_advect_thickness_x '//&
         '(1=no clip / smooth balanced, 0=full clip, 2=max compressive; faces where the limiter '//&
         'branch was not taken report 1.0)', 'nondim')
      CS%id_phi_y_FV = register_diag_field('ice_shelf_model','phi_y_FV',CS%diag%axesCv1, Time, &
         'Slope-limiter factor at each v-face from ice_shelf_advect_thickness_y '//&
         '(1=no clip / smooth balanced, 0=full clip, 2=max compressive; faces where the limiter '//&
         'branch was not taken report 1.0)', 'nondim')
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
      if (CS%use_DG_thickness .and. (CS%cism_taud == CISM_OFF)) then
        call calc_shelf_driving_stress_DG(CS, ISS, G, CS%taudx_shelf, CS%taudy_shelf)
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

  ! Without DG the nodal bed and thickness were only needed to derive the cell means.
  if (.not. CS%use_DG_thickness) then
    if (associated(CS%bed_node)) deallocate(CS%bed_node)
    if (associated(CS%h_nodal)) deallocate(CS%h_nodal)
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
                            ocean_mass, coupled_grounding, must_update_vel, vel_updated)
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
  logical,      optional, intent(out)   :: vel_updated !< True if the ice velocities were updated
                                              !! during this call.  This can not be anticipated by
                                              !! the caller, as the elapsed velocity time can also
                                              !! trigger an update.
  integer :: iters
  logical :: update_ice_vel, coupled_GL
  logical :: refresh_gfrac ! If true, refresh the grounded geometry after the advection, rather
                           ! than leaving it to the next velocity solve.

  update_ice_vel = .false.
  if (present(must_update_vel)) update_ice_vel = must_update_vel

  coupled_GL = .false.
  if (present(ocean_mass) .and. present(coupled_grounding)) coupled_GL = coupled_grounding

  ! The prescribed basal melt reads the grounded fraction, so it has to follow the ice thickness
  ! on every advective sub-step, not only on the sub-steps that solve for the velocities.  This
  ! is only ever true for the ice-only driver; ICE_ONLY_BASAL_MELT is not permitted otherwise.
  refresh_gfrac = (CS%ice_only_basal_melt .and. CS%advect_shelf) .and. (.not. coupled_GL)
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
  elseif (update_ice_vel .or. refresh_gfrac) then
    call update_OD_ffrac_uncoupled(CS, G, ISS%h_shelf(:,:))
    CS%GL_couple=.false.
  endif

  ! ice_shelf_solve_outer would otherwise do this itself, but the ice thickness does not change
  ! between here and there, so it is only ever done once per call.
  if (refresh_gfrac) then
    call update_grounded_geometry(CS, ISS, G)
    CS%grounded_geom_current = .true.
  endif

  if (update_ice_vel) then
    call ice_shelf_solve_outer(CS, ISS, G, US, CS%u_shelf, CS%v_shelf,CS%taudx_shelf,CS%taudy_shelf, iters, Time)
    CS%elapsed_velocity_time = 0.0
  endif

  if (present(vel_updated)) vel_updated = update_ice_vel

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
  ! DG(1) jump diagnostics. B-node quantities gather the co-located corners of the touching
  ! hmask=1 cells; face quantities compare the two sides at the face's endpoint nodes.
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_jump_n     ! Max-min corner thickness at a B-node [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_node_mx   ! Max corner thickness at a B-node [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_node_mn   ! Min corner thickness at a B-node [Z ~> m]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: h_jump_n_rel ! h_jump_n over the mean cell thickness [nondim]
  real :: h_cSW, h_cSE, h_cNW, h_cNE              ! Corner thicknesses at a B-node [Z ~> m]
  real :: Hb_SW, Hb_SE, Hb_NW, Hb_NE              ! Touching cells' mean thicknesses [Z ~> m]
  real :: hmax_b, hmin_b                          ! Max and min corner thickness [Z ~> m]
  real :: Hbar_sum                                ! Sum of the touching cell means [Z ~> m]
  real :: Hbar_avg                                ! Mean of the touching cell means [Z ~> m]
  integer :: n_valid                              ! Number of touching hmask=1 cells
  logical :: vSW, vSE, vNW, vNE                   ! True for touching hmask=1 cells
  integer :: ii, jj                               ! Touching-cell indices
  real, dimension(SZDIB_(G),SZDJ_(G)) :: hjump_fu ! u-face max|[h]| [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: sjump_fu ! u-face max|[s]| [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: sjump_fu_rel ! u-face max|[s]| over the face mean thickness [nondim]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: hjump_fu_sgn ! u-face mean signed [h] [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: sjump_fu_sgn ! u-face mean signed [s] [Z ~> m]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: un_fu    ! u-face |u| [L T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJB_(G)) :: hjump_fv ! v-face max|[h]| [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: sjump_fv ! v-face max|[s]| [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: sjump_fv_rel ! v-face max|[s]| over the face mean thickness [nondim]
  real, dimension(SZDI_(G),SZDJB_(G)) :: hjump_fv_sgn ! v-face mean signed [h] [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: sjump_fv_sgn ! v-face mean signed [s] [Z ~> m]
  real, dimension(SZDI_(G),SZDJB_(G)) :: un_fv    ! v-face |v| [L T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJ_(G)) :: eps_fu   ! u-face eps_e [T-1 ~> s-1]
  real, dimension(SZDI_(G),SZDJB_(G)) :: eps_fv   ! v-face eps_e [T-1 ~> s-1]
  real :: Hbar_face_avg                           ! Mean thickness of the two cells [Z ~> m]
  real :: rr                                      ! Ice to ocean density ratio [nondim]
  real :: bed1, bed2                              ! Bed at the face endpoints [Z ~> m]
  real :: h_m1, h_p1, h_m2, h_p2                  ! Minus/plus corner thickness at endpoints 1,2 [Z ~> m]
  real :: s_m1, s_p1, s_m2, s_p2                  ! Minus/plus corner surface at endpoints 1,2 [Z ~> m]

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
    if (CS%id_dg_damp_ahat > 0) call post_data(CS%id_dg_damp_ahat, CS%dg_damp_ahat, CS%diag)
    if (CS%id_dg_damp_spread > 0) call post_data(CS%id_dg_damp_spread, CS%dg_damp_spread, CS%diag)
    if (CS%id_dg_damp_tend > 0) call post_data(CS%id_dg_damp_tend, CS%dg_damp_tend, CS%diag)
    if (CS%id_dg_damp_gate > 0) call post_data(CS%id_dg_damp_gate, CS%dg_damp_gate, CS%diag)
    if (CS%id_dg_damp_want > 0) call post_data(CS%id_dg_damp_want, CS%dg_damp_want, CS%diag)
    if (CS%id_dg_damp_got > 0) call post_data(CS%id_dg_damp_got, CS%dg_damp_got, CS%diag)
    if (CS%id_dg_tilt_relax_rate > 0) &
      call post_data(CS%id_dg_tilt_relax_rate, CS%dg_tilt_relax_rate, CS%diag)
    if (CS%id_dg_tilt_relax_tend > 0) &
      call post_data(CS%id_dg_tilt_relax_tend, CS%dg_tilt_relax_tend, CS%diag)
    if (CS%id_dg_tilt_relax_alt_x > 0) &
      call post_data(CS%id_dg_tilt_relax_alt_x, CS%dg_tilt_relax_alt_x, CS%diag)
    if (CS%id_dg_tilt_relax_alt_y > 0) &
      call post_data(CS%id_dg_tilt_relax_alt_y, CS%dg_tilt_relax_alt_y, CS%diag)
    if (CS%id_dg_tilt_relax_alt_w > 0) &
      call post_data(CS%id_dg_tilt_relax_alt_w, CS%dg_tilt_relax_alt_w, CS%diag)
    if (CS%id_dg_tilt_relax_state > 0) &
      call post_data(CS%id_dg_tilt_relax_state, CS%dg_tilt_relax_state, CS%diag)
    if (CS%id_ground_frac > 0) call post_data(CS%id_ground_frac, CS%ground_frac, CS%diag)
    if (CS%id_f_ground_cell > 0) call post_data(CS%id_f_ground_cell, CS%f_ground_cell, CS%diag)
    if (CS%id_f_ground_node > 0) call post_data(CS%id_f_ground_node, CS%f_ground_node, CS%diag)
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
    if ((CS%id_h_jump_face_u > 0 .or. CS%id_h_jump_face_v > 0 .or. &
         CS%id_s_jump_face_u > 0 .or. CS%id_s_jump_face_v > 0 .or. &
         CS%id_s_jump_face_u_rel > 0 .or. CS%id_s_jump_face_v_rel > 0 .or. &
         CS%id_un_face_u > 0 .or. CS%id_un_face_v > 0) .and. associated(CS%h_nodal)) then
      ! Thickness and surface jumps across each face, with the face speed.
      call pass_corner_field(CS%h_nodal, G)
      call pass_var(CS%bed_node, G%domain, position=CORNER)
      call pass_vector(CS%u_shelf, CS%v_shelf, G%domain, TO_ALL, BGRID_NE)
      rr = CS%density_ice / CS%density_ocean_avg
      hjump_fu(:,:) = 0.0 ; sjump_fu(:,:) = 0.0 ; sjump_fu_rel(:,:) = 0.0 ; un_fu(:,:) = 0.0
      hjump_fu_sgn(:,:) = 0.0 ; sjump_fu_sgn(:,:) = 0.0
      hjump_fv(:,:) = 0.0 ; sjump_fv(:,:) = 0.0 ; sjump_fv_rel(:,:) = 0.0 ; un_fv(:,:) = 0.0
      hjump_fv_sgn(:,:) = 0.0 ; sjump_fv_sgn(:,:) = 0.0
      ! u-faces: minus = west cell, plus = east cell; node 1 south, 2 north.
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
      ! v-faces: minus = south cell, plus = north cell; node 1 west, 2 east.
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
    ! Face eps_e with the art-visc stencil, posted even when art visc is off.
    if (CS%id_dg_eps_face_u > 0) then
      eps_fu(:,:) = 0.0
      do j = G%jsc, G%jec ; do I = G%IscB, G%IecB
        if (ISS%hmask(I,  j) /= 1.0 .and. ISS%hmask(I,  j) /= 3.0) cycle
        if (ISS%hmask(I+1,j) /= 1.0 .and. ISS%hmask(I+1,j) /= 3.0) cycle
        eps_fu(I,j) = dg1_eps_face_u(CS, G, I, j)
      enddo ; enddo
      call post_data(CS%id_dg_eps_face_u, eps_fu, CS%diag)
    endif
    if (CS%id_dg_eps_face_v > 0) then
      eps_fv(:,:) = 0.0
      do J = G%JscB, G%JecB ; do i = G%isc, G%iec
        if (ISS%hmask(i,J  ) /= 1.0 .and. ISS%hmask(i,J  ) /= 3.0) cycle
        if (ISS%hmask(i,J+1) /= 1.0 .and. ISS%hmask(i,J+1) /= 3.0) cycle
        eps_fv(i,J) = dg1_eps_face_v(CS, G, i, J)
      enddo ; enddo
      call post_data(CS%id_dg_eps_face_v, eps_fv, CS%diag)
    endif
    if (CS%id_h_source_rate > 0 .and. associated(CS%h_source_rate_last)) &
        call post_data(CS%id_h_source_rate, CS%h_source_rate_last, CS%diag)
    if (CS%id_phi_x_FV > 0 .and. associated(CS%phi_x_FV)) &
        call post_data(CS%id_phi_x_FV, CS%phi_x_FV, CS%diag)
    if (CS%id_phi_y_FV > 0 .and. associated(CS%phi_y_FV)) &
        call post_data(CS%id_phi_y_FV, CS%phi_y_FV, CS%diag)
    if (CS%id_dg_art_visc_coef_u > 0 .and. associated(CS%dg_art_visc_coef_u)) &
        call post_data(CS%id_dg_art_visc_coef_u, CS%dg_art_visc_coef_u, CS%diag)
    if (CS%id_dg_art_visc_coef_v > 0 .and. associated(CS%dg_art_visc_coef_v)) &
        call post_data(CS%id_dg_art_visc_coef_v, CS%dg_art_visc_coef_v, CS%diag)
    if (CS%id_dg_art_visc_nu_u > 0 .and. associated(CS%dg_art_visc_nu_u)) &
        call post_data(CS%id_dg_art_visc_nu_u, CS%dg_art_visc_nu_u, CS%diag)
    if (CS%id_dg_art_visc_nu_v > 0 .and. associated(CS%dg_art_visc_nu_v)) &
        call post_data(CS%id_dg_art_visc_nu_v, CS%dg_art_visc_nu_v, CS%diag)
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

  ! Hold the prescribed thickness on boundary cells.  Under DG the cell carries its own corner
  ! values, and recompute_h_shelf_from_nodal already publishes their mean, so there is nothing to
  ! impose here; h_bdry_val is the fallback for the non-nodal path.  The test is on the mask, which
  ! is what the flux, the viscosity and the front term all use.
  do j=jsd,jed ; do i=isd,ied ; if (ISS%hmask(i,j) == 3.0) then
    if (.not.CS%use_DG_thickness) ISS%h_shelf(i,j) = CS%h_bdry_val(i,j)
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

!> Refresh every grounded/floating geometry field that is a function of the current ice thickness:
!! CS%H_node, the flotation gate and corner flotation fields, CS%ground_frac (along with
!! CS%xi_basal), and the analytic quadrant grounded fractions
!! CS%f_ground_cell and CS%f_ground_node.  These are read by the driving stress and the basal
!! friction in the velocity solve, and by the prescribed basal melt parameterization.
subroutine update_grounded_geometry(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe
                                               !! the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.

  real :: rhoi_rhow ! The density of ice divided by a typical water density [nondim]
  integer :: i, j

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  ! need to make these conditional on GL interpolation
  CS%H_node(:,:) = 0.0
  !CS%ground_frac(:,:) = 0.0

  if (.not. CS%GL_couple) then
    do j=G%jsc,G%jec ; do i=G%isc,G%iec
      if (rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf) - CS%bed_elev(i,j) > 0) then
        CS%ground_frac(i,j) = 1.0
        CS%OD_av(i,j) =0.0
      endif
    enddo ; enddo
  endif

  ! Fractional ground_frac in GL_regularize cells, before the driving stress uses it.
  ! H_node is read by the non-DG compute_ground_frac and CG_action_subgrid_basal.
  if (CS%GL_regularize .and. .not. CS%use_DG_thickness) then
    call interpolate_H_to_B(G, ISS%h_shelf, ISS%hmask, CS%H_node, CS%min_h_shelf)
  endif
  ! FV sub-element corner fields, so ground_frac uses the same flotation field.
  if (CS%i2n_friction .or. CS%i2n_taud) &
    call build_corner_flotation_fields(CS, ISS, G)
  call compute_ground_frac(CS, ISS, G, CS%H_node)

  ! Analytic quadrant grounding-line fractions for friction and/or the driving-stress surface
  ! blend (Leguy et al. 2021). Uses cell-mean h_shelf/bed_elev, so it is independent of the
  ! thickness-advection scheme.
  if (CS%gl_quad_friction) call compute_gl_quadrant_fractions(CS, ISS, G)

end subroutine update_grounded_geometry

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
  real, dimension(SZDIB_(G),SZDJB_(G)) :: u_pre_newton, v_pre_newton ! Velocities saved at the
                                              ! Picard-to-Newton switch, restored if Newton
                                              ! diverges [L T-1 ~> m s-1]
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
  real    :: newton_after_tol_loc ! Working Picard-to-Newton switch threshold for this solve
                                  ! [nondim]; reduced tenfold on each divergence rescue.
  real    :: err_newton_enter ! Nonlinear residual at the Picard-to-Newton switch, the
                              ! reference for the divergence test (same units as err_max)
  real    :: Norm_newton_enter ! Norm saved at the switch for restore (err mode 3) [L T-1 ~> m s-1]
  logical :: rescue_enabled   ! Newton divergence rescue is configured and applicable
  logical :: newton_armed     ! Pre-Newton state is saved; divergence rescue available
  logical :: diverging        ! The current Newton residual triggers a rescue
  integer :: n_rescue         ! Number of divergence rescues performed in this solve

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

  ! Warning: This turns off Picard entirely and may not converge.
  if (CS%newton_after_tolerance<=0.0) CS%doing_newton=.true.

  ! Refresh the grounded geometry that the driving stress and basal friction integrate over.
  ! h_shelf has not changed since update_ice_shelf did this, in the cases where it does, so
  ! the work is not repeated here.
  if (.not. CS%grounded_geom_current) call update_grounded_geometry(CS, ISS, G)
  CS%grounded_geom_current = .false.

  ! Calculate RHS. CISM_TAUD uses the FV driving stress on the cell means, even under DG.
  if (CS%use_DG_thickness .and. (CS%cism_taud == CISM_OFF)) then
    call calc_shelf_driving_stress_DG(CS, ISS, G, taudx, taudy)
  else
    call calc_shelf_driving_stress(CS, ISS, G, US, taudx, taudy, CS%OD_av)
  endif
  call pass_vector(taudx, taudy, G%domain, TO_ALL, BGRID_NE)

  ! Calculate basal drag constants and initial velocity
  call calc_shelf_basal_prefactors(CS, ISS, G, US)
  if (CS%local_basal_friction) call calc_shelf_basal_prefactors_node(CS, ISS, G, US)
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
    call CG_action(CS, Au, Av, u_shlf, v_shlf, CS%Phi, CS%Phisub, CS%umask, CS%vmask, ISS%hmask, CS%H_node, &
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

  ! Newton divergence rescue state. The working switch threshold is per-solve: rescues
  ! reduce it tenfold, and the configured value is restored at the next solve.
  rescue_enabled = CS%newton_divergence_rescue .and. (CS%newton_after_tolerance > 0.0)
  newton_after_tol_loc = CS%newton_after_tolerance
  newton_armed = .false.
  err_newton_enter = 0.0 ; Norm_newton_enter = 0.0
  n_rescue = 0

  !! begin loop

  iter = 0
  do
    iter = iter + 1
    if (iter > 50) exit

    ! The linear solve
    call ice_shelf_solve_inner(CS, ISS, G, US, u_shlf, v_shlf, taudx, taudy, CS%H_node, &
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
        CS%H_node, CS%ice_visc, CS%bed_elev, u_shlf, v_shlf, &
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

      ! Newton divergence rescue: if the residual has gone NaN or grown beyond
      ! NEWTON_DIVERGENCE_FACTOR times its value at the Picard-to-Newton switch,
      ! Newton was activated outside its basin of attraction. Restore the
      ! pre-Newton iterate, revert to Picard with a fresh outer-iteration budget,
      ! and lower the working switch threshold tenfold. Rescues self-terminate:
      ! once the threshold drops below the convergence tolerance Newton cannot
      ! re-activate, so the solve completes as pure Picard.
      if (rescue_enabled .and. CS%doing_newton .and. newton_armed) then
        diverging = (err_max /= err_max) .or. &
                    (err_max > CS%newton_divergence_factor * err_newton_enter)
        if (diverging) then
          n_rescue = n_rescue + 1
          u_shlf(:,:) = u_pre_newton(:,:) ; v_shlf(:,:) = v_pre_newton(:,:)
          u_last(:,:) = u_shlf(:,:) ; v_last(:,:) = v_shlf(:,:)
          if (CS%nonlin_solve_err_mode == 3) Norm = Norm_newton_enter
          CS%doing_newton = .false. ; newton_armed = .false.
          ew_prev_resid = 0.0
          CS%cg_tol_current = CS%cg_tolerance
          if (n_rescue >= CS%newton_max_rescues) then
            ! Rescue budget exhausted: pure Picard for the remainder of this solve.
            newton_after_tol_loc = 0.0
          else
            newton_after_tol_loc = 0.1 * newton_after_tol_loc
          endif
          iter = 0
          ! Rebuild the (Picard) viscosity of the restored iterate so the next
          ! inner solve does not reuse operators from the divergent state.
          call calc_shelf_visc(CS, ISS, G, US, u_shlf, v_shlf)
          call pass_var(CS%ice_visc, G%domain, complete=.true.)
          write(mesg,*) "ice_shelf_solve_outer: Newton diverged (rescue ", n_rescue, &
              "); restored pre-Newton state, switch threshold now ", newton_after_tol_loc
          call MOM_mesg(mesg, 2)
          cycle
        endif
      endif

      ! Activate Newton
      if (err_max <= newton_after_tol_loc * err_init .and. .not. CS%doing_newton) then
        if (rescue_enabled) then
          ! Save the switch state so a divergent Newton excursion can be undone.
          u_pre_newton(:,:) = u_shlf(:,:) ; v_pre_newton(:,:) = v_shlf(:,:)
          err_newton_enter = err_max
          if (CS%nonlin_solve_err_mode == 3) Norm_newton_enter = Norm
          newton_armed = .true.
        endif
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
              CS%H_node, CS%ice_visc, CS%bed_elev, u_shlf, v_shlf, &
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
  real :: cfl_wt     ! The Lax-Wendroff (1-CFL) reconstruction weight, or 1 if disabled [nondim]

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
            slope_lim = slope_limiter(h0(i,j)-h0(i-1,j), h0(i+1,j)-h0(i,j), CS%adv_thickness_limiter)
            ! This is a 2nd-order scheme with a TVD slope limiter, optionally Lax-Wendroff (1-CFL)
            ! weighted for time accuracy.  We could try PPM here.
            cfl_wt = 1.0
            if (CS%adv_cfl_weight) cfl_wt = max(0.0, 1.0 - u_face*time_step*G%IdxT(i,j))
            h_face = h0(i,j) - slope_lim * (0.5 * cfl_wt * (h0(i,j)-h0(i+1,j)))
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
            slope_lim = slope_limiter(h0(i+1,j)-h0(i,j), h0(i+2,j)-h0(i+1,j), CS%adv_thickness_limiter)
            cfl_wt = 1.0
            if (CS%adv_cfl_weight) cfl_wt = max(0.0, 1.0 + u_face*time_step*G%IdxT(i+1,j))
            h_face = h0(i+1,j) - slope_lim * (0.5 * cfl_wt * (h0(i+2,j)-h0(i+1,j)))
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
  real :: cfl_wt     ! The Lax-Wendroff (1-CFL) reconstruction weight, or 1 if disabled [nondim]

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
            slope_lim = slope_limiter(h0(i,j)-h0(i,j-1), h0(i,j+1)-h0(i,j), CS%adv_thickness_limiter)
            ! This is a 2nd-order scheme with a TVD slope limiter, optionally Lax-Wendroff (1-CFL)
            ! weighted for time accuracy.  We could try PPM here.
            cfl_wt = 1.0
            if (CS%adv_cfl_weight) cfl_wt = max(0.0, 1.0 - v_face*time_step*G%IdyT(i,j))
            h_face = h0(i,j) - slope_lim * (0.5 * cfl_wt * (h0(i,j)-h0(i,j+1)))
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
            slope_lim = slope_limiter(h0(i,j+1)-h0(i,j), h0(i,j+2)-h0(i,j+1), CS%adv_thickness_limiter)
            cfl_wt = 1.0
            if (CS%adv_cfl_weight) cfl_wt = max(0.0, 1.0 + v_face*time_step*G%IdyT(i,j+1))
            h_face = h0(i,j+1) - slope_lim * (0.5 * cfl_wt * (h0(i,j+2)-h0(i,j+1)))
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
              ! Keep the DG corners equal to the new cell mean.
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
  logical, dimension(SZDI_(G),SZDJ_(G)) :: grnd ! True at grounded ice cells, set by the same
                       ! flotation test that builds S; used by the one-sided grounding-line driving stress.
  real    :: rho, rhow, rhoi_rhow ! Ice and ocean densities [R ~> kg m-3]
  real    :: sx, sy    ! Ice shelf top slopes at tracer points [Z L-1 ~> nondim]
  real    :: hx, hy    ! Effective ice thickness for the one-sided grounding-line driving stress
                       ! in the x- and y- directions, averaged across the face used [Z ~> m].
  logical :: gnd_E, gnd_W, gnd_N, gnd_S ! True if the neighboring cell is grounded ice.
  logical :: flt_E, flt_W, flt_N, flt_S ! True if the neighboring cell is floating ice.
  logical :: gnd_EE, gnd_WW, gnd_NN, gnd_SS ! As above, but for the next cell out.
  logical :: flt_EE, flt_WW, flt_NN, flt_SS ! As above, but for the next cell out.
  logical :: this_gnd, this_flt ! True if the current cell is grounded / floating ice.
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


  ! Sub-element driving stress on the same grounding-line partition the friction integrates over;
  ! hand off and return so every caller routes through it.
  if (CS%i2n_taud) then
    call calc_shelf_driving_stress_i2n(CS, ISS, G, US, taudx, taudy)
    return
  endif

  ! Compact nodal surface-gradient driving stress (Lipscomb et al. 2019 eq. 14) is a drop-in
  ! alternative for the FV path; hand off and return so every caller routes through it.
  if (CS%fv_taud_vertex_grad) then
    call calc_shelf_driving_stress_vertex(CS, ISS, G, US, taudx, taudy, OD)
    return
  endif

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

  ! Smooth the surface across the grounding line using the analytic cell grounded fraction
  ! (Leguy et al. 2021). Mutually exclusive with FV_GL_ONE_SIDED_TAUD (enforced at init).

  call pass_var(S, G%domain)

  ! Flag grounded ice cells for the one-sided grounding-line driving stress, using the same
  ! cell-center flotation test that builds S (Cornford et al. 2013, Feldmann et al. 2014). The
  ! mask is built over the full data domain so the +/-2-cell grounding-line stencil can be read
  ! directly; bed_elev, h_shelf, hmask and ground_frac already carry valid halo values here.
  if (CS%FV_GL_one_sided) then
    do j=jsd,jed ; do i=isd,ied
      if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then
        if (CS%GL_couple) then
          grnd(i,j) = (CS%ground_frac(i,j) >= 1.0)
        else
          grnd(i,j) = (rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf) - CS%bed_elev(i,j) > 0.0)
        endif
      else
        !Non-ice cells for the grounding-line stencil are treated as grounded/floating cells
        !according to whether they are land/ocean cells
        if (CS%bed_elev(i,j)>0) then
          grnd(i,j) = .false.
        else
          grnd(i,j) = .true.
        endif
      endif
    enddo ; enddo
  endif

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

        ! Effective thickness for the driving stress; replaced by a face-averaged value below
        ! wherever the one-sided grounding-line treatment is applied.
        hx = max(ISS%h_shelf(i,j),CS%min_h_shelf)
        hy = max(ISS%h_shelf(i,j),CS%min_h_shelf)

        ! One-sided finite-volume driving stress near the grounding line, following
        ! Cornford et al. (2013) eqs 27-29 (and Feldmann et al. 2014). Wherever a centered
        ! slope would straddle the grounding line, it is replaced by a one-sided difference
        ! that stays on a single side (grounded or floating) of the line, with the thickness
        ! averaged across the face used. Grounded/floating is the cell-center flotation test
        ! held in grnd
        if (CS%FV_GL_one_sided) then
          ! x-direction. Require the +/-2 stencil to lie within the data domain.
          if ((i-2 >= isd) .and. (i+2 <= ied)) then
            this_gnd = grnd(i,j)
            this_flt = .not. this_gnd
            gnd_E  = grnd(i+1,j) ; flt_E = valid_E .and. (.not. gnd_E)
            gnd_W  = grnd(i-1,j) ; flt_W = valid_W .and. (.not. gnd_W)
            gnd_EE = grnd(i+2,j) ; flt_EE = (.not. gnd_EE)
            gnd_WW = grnd(i-2,j) ; flt_WW = (.not. gnd_WW)

            if (this_gnd .and. gnd_W .and. flt_E .and. flt_EE) then
              ! Last grounded cell, grounding line to the east (eqs 27-28): backward difference.
              sx = (S(i,j) - S(i-1,j)) * G%IdxCu(I-1,j)
              hx = 0.5*(max(ISS%h_shelf(i,j),CS%min_h_shelf) + max(ISS%h_shelf(i-1,j),CS%min_h_shelf))
            elseif (this_gnd .and. gnd_E .and. flt_W .and. flt_WW) then
              ! Last grounded cell, grounding line to the west: forward difference.
              sx = (S(i+1,j) - S(i,j)) * G%IdxCu(I,j)
              hx = 0.5*(max(ISS%h_shelf(i+1,j),CS%min_h_shelf) + max(ISS%h_shelf(i,j),CS%min_h_shelf))
            elseif (this_flt .and. flt_E .and. gnd_W .and. gnd_WW) then
              ! First floating cell, grounding line to the west (eq 29): forward difference.
              sx = (S(i+1,j) - S(i,j)) * G%IdxCu(I,j)
              hx = 0.5*(max(ISS%h_shelf(i+1,j),CS%min_h_shelf) + max(ISS%h_shelf(i,j),CS%min_h_shelf))
            elseif (this_flt .and. flt_W .and. gnd_E .and. gnd_EE) then
              ! First floating cell, grounding line to the east: backward difference.
              sx = (S(i,j) - S(i-1,j)) * G%IdxCu(I-1,j)
              hx = 0.5*(max(ISS%h_shelf(i,j),CS%min_h_shelf) + max(ISS%h_shelf(i-1,j),CS%min_h_shelf))
            endif
          endif

          ! y-direction. Require the +/-2 stencil to lie within the data domain.
          if ((j-2 >= jsd) .and. (j+2 <= jed)) then
            this_gnd = grnd(i,j)
            this_flt = .not. this_gnd
            gnd_N  = grnd(i,j+1) ; flt_N = valid_N .and. (.not. gnd_N)
            gnd_S  = grnd(i,j-1) ; flt_S = valid_S .and. (.not. gnd_S)
            gnd_NN = grnd(i,j+2) ; flt_NN = (.not. gnd_NN)
            gnd_SS = grnd(i,j-2) ; flt_SS = (.not. gnd_SS)

            if (this_gnd .and. gnd_S .and. flt_N .and. flt_NN) then
              ! Last grounded cell, grounding line to the north (eqs 27-28): backward difference.
              sy = (S(i,j) - S(i,j-1)) * G%IdyCv(i,J-1)
              hy = 0.5*(max(ISS%h_shelf(i,j),CS%min_h_shelf) + max(ISS%h_shelf(i,j-1),CS%min_h_shelf))
            elseif (this_gnd .and. gnd_N .and. flt_S .and. flt_SS) then
              ! Last grounded cell, grounding line to the south: forward difference.
              sy = (S(i,j+1) - S(i,j)) * G%IdyCv(i,J)
              hy = 0.5*(max(ISS%h_shelf(i,j+1),CS%min_h_shelf) + max(ISS%h_shelf(i,j),CS%min_h_shelf))
            elseif (this_flt .and. flt_N .and. gnd_S .and. gnd_SS) then
              ! First floating cell, grounding line to the south (eq 29): forward difference.
              sy = (S(i,j+1) - S(i,j)) * G%IdyCv(i,J)
              hy = 0.5*(max(ISS%h_shelf(i,j+1),CS%min_h_shelf) + max(ISS%h_shelf(i,j),CS%min_h_shelf))
            elseif (this_flt .and. flt_S .and. gnd_N .and. gnd_NN) then
              ! First floating cell, grounding line to the north: backward difference.
              sy = (S(i,j) - S(i,j-1)) * G%IdyCv(i,J-1)
              hy = 0.5*(max(ISS%h_shelf(i,j),CS%min_h_shelf) + max(ISS%h_shelf(i,j-1),CS%min_h_shelf))
            endif
          endif
        endif

        if (CS%max_surface_slope>0) then
          scale = CS%max_surface_slope / max( sqrt((sx**2) + (sy**2)), CS%max_surface_slope )
          sx = scale*sx ; sy = scale*sy
        endif

        sx_e(i,j) = (-.25 * G%areaT(i,j)) * ((rho * grav) * (hx * sx))
        sy_e(i,j) = (-.25 * G%areaT(i,j)) * ((rho * grav) * (hy * sy))

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

!> Finite-volume (non-DG) driving stress with the surface gradient evaluated directly at B-grid
!! nodes from the four surrounding cell-center surface elevations, following Lipscomb et al. (2019,
!! Geosci. Model Dev. 12:387-424) eq. 14 and their "option 3" ice-margin treatment. This compact
!! 4-cell stencil replaces the wider cell-centroid centered difference used by
!! calc_shelf_driving_stress, and so is less smeared across the grounding line. Per parallel edge,
!! a gradient is included when both cells are ice-covered, or when an ice-covered cell lies above an
!! ice-free land neighbor; edges across an ice/ice-free-ocean margin contribute no gradient (the
!! lateral pressure there is supplied by the Neumann face term, retained verbatim below), and
!! nunatak edges (ice below ice-free land) likewise contribute nothing. The body force is then
!! assembled by element quadrature (A4, eqs A26-A27): the nodal slope is interpolated to each 2x2
!! Gauss point with the bilinear basis and distributed to the corner nodes weighted by the basis,
!! the same consistent integration the basal friction uses in CG_action, so the driving-stress and
!! friction grounding lines co-locate. The surface field S is built exactly as in
!! calc_shelf_driving_stress, so this honors
!! MAX_SURFACE_SLOPE. Selected by CISM_TAUD; mutually exclusive with
!! FV_GL_ONE_SIDED_TAUD (a centroid-slope construct with no nodal analog).
subroutine calc_shelf_driving_stress_vertex(CS, ISS, G, US, taudx, taudy, OD)
  type(ice_shelf_dyn_CS), intent(in)   :: CS  !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(in)    :: ISS !< A structure describing the ice-shelf state
  type(ocean_grid_type), intent(inout) :: G   !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)    :: US  !< A structure containing unit conversion factors
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in)    :: OD  !< ocean floor depth at tracer points [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudx  !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudy  !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]

  real, dimension(SZDI_(G),SZDJ_(G))   :: S     ! surface elevation [Z ~> m].
  real, dimension(SZDIB_(G),SZDJB_(G)) :: sx_n, sy_n ! Nodal surface slopes [Z L-1 ~> nondim]
  logical, dimension(SZDI_(G),SZDJ_(G)) :: ice_cell  ! True at ice-covered cells (grounded or floating)
  logical, dimension(SZDI_(G),SZDJ_(G)) :: land_cell ! True at ice-free cells at/above sea level (bed_elev<=0)
  real    :: rho, rhow, rhoi_rhow ! Ice and ocean densities [R ~> kg m-3] and their ratio [nondim]
  real    :: grav      ! The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real    :: num_x, num_y ! Sum of the valid parallel-edge surface differences at a node [Z ~> m]
  real    :: den_x, den_y ! Sum of the valid parallel-edge lengths at a node [L ~> m]
  real    :: g_x1, g_x2   ! Per-edge x-direction surface slopes at a node; 0 on a dropped edge [Z L-1 ~> nondim]
  real    :: g_y1, g_y2   ! Per-edge y-direction surface slopes at a node; 0 on a dropped edge [Z L-1 ~> nondim]
  integer :: ndom_x, ndom_y ! Number of in-domain parallel edges at a node, the gradient divisor
  integer :: nok_x, nok_y   ! Number of those edges that also carry ice and so contribute
  real    :: smag      ! Surface slope magnitude at a node [Z L-1 ~> nondim]
  real    :: scale     ! Scaling factor enforcing MAX_SURFACE_SLOPE [nondim]
  real    :: neumann_val ! Lateral-pressure boundary term [R Z L2 T-2 ~> kg s-2]
  real    :: xquad(2)  ! 2-point Gauss-Legendre quadrature locations on [0,1] [nondim]
  real    :: He        ! Cell-mean ice thickness used in the driving-stress integral [Z ~> m]
  real    :: wq        ! Per-quadrature-point area weight, 1/4 |J_q| (= 1/4 areaT for rectangular cells) [L2 ~> m2]
  real    :: dsdx_qp, dsdy_qp ! Surface slope interpolated to a quadrature point [Z L-1 ~> nondim]
  real    :: phim(2,2) ! Bilinear nodal basis values of the 4 cell corners at a quadrature point [nondim]
  real    :: hA        ! h-weighted node control mass for the lumped assembly, sum of the metric
                       ! lumped corner masses (1/4 |J_q|-integrated basis) over the ice-covered
                       ! cells around the node [Z L2 ~> m3]
  real    :: hA_sw, hA_se, hA_nw, hA_ne ! Per-cell metric-lumped corner-mass weights at the four cells
                       ! around a node, summed in diagonal pairs for rotation invariance [Z L2 ~> m3]
  real    :: txqp(2,2,4), tyqp(2,2,4) ! Per-quadrature-point nodal driving-stress contributions within
                       ! one element, pair-summed over quadrature points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: taudx_b, taudy_b ! Per-element corner contributions to the
                       ! nodal driving stress, diagonal-pair-summed for rotation invariance [R L3 Z T-2 ~> kg m s-2]
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed
  integer :: iq, jq, iphi, jphi, ilq, jlq, Itgt, Jtgt, qp
  integer :: i_off, j_off, gisc, gjsc, giec, gjec

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset
  gisc = 1 ; gjsc = 1 ; giec = G%domain%niglobal ; gjec = G%domain%njglobal

  rho = CS%density_ice ; rhow = CS%density_ocean_avg ; grav = CS%g_Earth
  rhoi_rhow = rho/rhow
  xquad(1) = .5*(1. - sqrt(1./3.)) ; xquad(2) = .5*(1. + sqrt(1./3.))

  ! Surface elevation S -- identical to calc_shelf_driving_stress.
  if (CS%GL_couple) then
    do j=jsc-2,jec+2 ; do i=isc-2,iec+2
      S(i,j) = -CS%bed_elev(i,j) + (OD(i,j) + max(ISS%h_shelf(i,j),CS%min_h_shelf))
    enddo ; enddo
  else
    do j=jsc-2,jec+2 ; do i=isc-2,iec+2
      if (rhoi_rhow * max(ISS%h_shelf(i,j),CS%min_h_shelf) - CS%bed_elev(i,j) <= 0) then
        S(i,j) = (1 - rhoi_rhow)*max(ISS%h_shelf(i,j),CS%min_h_shelf)
      else
        S(i,j) = max(ISS%h_shelf(i,j),CS%min_h_shelf)-CS%bed_elev(i,j)
      endif
    enddo ; enddo
  endif
  call pass_var(S, G%domain)

  ! Cell classification for the option-3 margin rule: ice-covered (hmask 1 or 3), ice-free land
  ! (no ice and bed at/above sea level, bed_elev<=0), else ice-free ocean.
  do j=jsd,jed ; do i=isd,ied
    ice_cell(i,j)  = (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3)
    land_cell(i,j) = (.not. ice_cell(i,j)) .and. (CS%bed_elev(i,j) <= 0.0)
  enddo ; enddo

  ! Nodal surface gradient (generalized Lipscomb 2019 eq. 14), built per parallel edge so margin
  ! edges can be dropped (option 3). Node (I,J) is the NE corner of cell (i,j); its four cells are
  ! (i,j),(i+1,j),(i,j+1),(i+1,j+1) (uppercase I,J equal lowercase i,j in Fortran). Computed over a
  ! node range wide enough to cover every corner of the element-integration loop below.
  !
  ! When both parallel edges are usable, they are combined as a ratio of summed surface differences
  ! over summed edge lengths -- num/den -- which is the isoparametric bilinear gradient CG_action
  ! forms for nodal fields (metric-interpolated numerator over metric-interpolated denominator), not
  ! an average of the two per-edge ratios. The two are bitwise identical on a uniform grid but differ
  ! where the parallel edges have unequal length (e.g. lat/lon, where dxCu varies with latitude), so
  ! this matches CG_action's metric treatment.
  !
  ! When an edge is dropped, the two reasons are treated differently, following CISM
  ! (glissade_surface_elevation_gradient with ho_gradient = HO_GRADIENT_CENTERED, which forms
  ! 0.5*(ds_dx_edge(i,j) + ds_dx_edge(i,j+1)) with a masked edge left at zero):
  !
  ! (a) The neighbor cell is genuinely ice-free (the option-3 margin rule). CISM leaves that edge's
  !     gradient at zero but still divides by two, so the nodal slope is HALVED. The surviving edge
  !     therefore carries a fixed 1/2 weight here, not a renormalized full weight.
  ! (b) The edge crosses a non-reentrant domain wall. CISM has no such case: its halo BC is periodic
  !     (parallel_halo), so the across-wall cells hold ice and both edges are valid. The faithful
  !     analogue is a mirror BC, under which the out-of-domain edge gradient equals the in-domain one
  !     and 0.5*(g+g) = g -- i.e. renormalize over the in-domain edges. Halving here instead would
  !     cut the wall-node driving stress in half while its friction and viscous terms keep their
  !     (half) control volume, breaking the meridional symmetry of channel configs (MISMIP3D).
  !
  ! So the divisor is the number of IN-DOMAIN parallel edges, while only the ice-usable ones
  ! contribute. With no ice-free cells anywhere (e.g. MISMIP3D, fixed front, hmask=1 everywhere)
  ! case (a) never arises and this reduces exactly to the renormalized num/den form.
  sx_n(:,:) = 0.0 ; sy_n(:,:) = 0.0
  do j=jsc-2,jec+1 ; do i=isc-2,iec+1
    ! x-slope: south edge (i,j)->(i+1,j) and north edge (i,j+1)->(i+1,j+1)
    num_x = 0.0 ; den_x = 0.0 ; g_x1 = 0.0 ; g_x2 = 0.0 ; ndom_x = 0 ; nok_x = 0
    if (edge_in_domain(i,j,  i+1,j  )) then
      ndom_x = ndom_x + 1
      if (edge_ice_ok(i,j,  i+1,j  )) then
        nok_x = nok_x + 1 ; g_x1 = (S(i+1,j)   - S(i,j)  ) / G%dxCu(I,j)
        num_x = num_x + (S(i+1,j)   - S(i,j)  ) ; den_x = den_x + G%dxCu(I,j)
      endif
    endif
    if (edge_in_domain(i,j+1,i+1,j+1)) then
      ndom_x = ndom_x + 1
      if (edge_ice_ok(i,j+1,i+1,j+1)) then
        nok_x = nok_x + 1 ; g_x2 = (S(i+1,j+1) - S(i,j+1)) / G%dxCu(I,j+1)
        num_x = num_x + (S(i+1,j+1) - S(i,j+1)) ; den_x = den_x + G%dxCu(I,j+1)
      endif
    endif
    if (nok_x == ndom_x) then
      if (den_x > 0.0) sx_n(I,J) = num_x / den_x
    elseif (nok_x > 0) then
      sx_n(I,J) = (g_x1 + g_x2) / real(ndom_x)
    endif

    ! y-slope: west edge (i,j)->(i,j+1) and east edge (i+1,j)->(i+1,j+1)
    num_y = 0.0 ; den_y = 0.0 ; g_y1 = 0.0 ; g_y2 = 0.0 ; ndom_y = 0 ; nok_y = 0
    if (edge_in_domain(i,j,  i,  j+1)) then
      ndom_y = ndom_y + 1
      if (edge_ice_ok(i,j,  i,  j+1)) then
        nok_y = nok_y + 1 ; g_y1 = (S(i,j+1)   - S(i,j)  ) / G%dyCv(i,J)
        num_y = num_y + (S(i,j+1)   - S(i,j)  ) ; den_y = den_y + G%dyCv(i,J)
      endif
    endif
    if (edge_in_domain(i+1,j,i+1,j+1)) then
      ndom_y = ndom_y + 1
      if (edge_ice_ok(i+1,j,i+1,j+1)) then
        nok_y = nok_y + 1 ; g_y2 = (S(i+1,j+1) - S(i+1,j)) / G%dyCv(i+1,J)
        num_y = num_y + (S(i+1,j+1) - S(i+1,j)) ; den_y = den_y + G%dyCv(i+1,J)
      endif
    endif
    if (nok_y == ndom_y) then
      if (den_y > 0.0) sy_n(I,J) = num_y / den_y
    elseif (nok_y > 0) then
      sy_n(I,J) = (g_y1 + g_y2) / real(ndom_y)
    endif

    ! Cap the surface slope magnitude (MAX_SURFACE_SLOPE), as in calc_shelf_driving_stress.
    if (CS%max_surface_slope > 0) then
      smag = sqrt((sx_n(I,J)**2) + (sy_n(I,J)**2))
      scale = CS%max_surface_slope / max(smag, CS%max_surface_slope)
      sx_n(I,J) = scale*sx_n(I,J) ; sy_n(I,J) = scale*sy_n(I,J)
    endif
  enddo ; enddo

  if (CS%local_fv_taud_vertex) then
    ! Mass-lumped ("local") assembly: each node's driving stress uses its own slope alone,
    ! tau_d = -rho g grad(s) * hA, with hA the h-weighted lumped nodal mass over the ice-covered
    ! cells around the node (Lipscomb 2019 A4 local method; more robust for sharp surfaces).
    do j=jsc-1,jec ; do i=isc-1,iec
      ! Metric-exact lumped (row-sum) mass contribution of each of the four cells around node (I,J):
      ! h * integral of that cell's node-corner bilinear basis against the element Jacobian,
      ! (h/4) sum_qp phi_corner(qp) Jac(qp) (lumped_corner_mass below; 2-pt Gauss is exact for the
      ! bilinear-orthogonal map, so this carries the lat/lon metric exactly and reduces to
      ! 0.25*areaT*h on rectangular cells). The corner passed is the cell's corner that coincides
      ! with the node: cell (i,j) -> NE, (i+1,j) -> NW, (i,j+1) -> SE, (i+1,j+1) -> SW. Cells are
      ! summed in diagonal pairs (SW+NE)+(SE+NW) so the lumped mass is bitwise invariant under
      ! 90 deg rotation, matching CG_action.
      hA_sw = 0.0 ; if (ice_cell(i,  j  )) hA_sw = lumped_corner_mass(i,  j,   2, 2)
      hA_se = 0.0 ; if (ice_cell(i+1,j  )) hA_se = lumped_corner_mass(i+1,j,   1, 2)
      hA_nw = 0.0 ; if (ice_cell(i,  j+1)) hA_nw = lumped_corner_mass(i,  j+1, 2, 1)
      hA_ne = 0.0 ; if (ice_cell(i+1,j+1)) hA_ne = lumped_corner_mass(i+1,j+1, 1, 1)
      hA = (hA_sw + hA_ne) + (hA_se + hA_nw)
      taudx(I,J) = taudx(I,J) - (rho*grav) * (hA * sx_n(I,J))
      taudy(I,J) = taudy(I,J) - (rho*grav) * (hA * sy_n(I,J))
    enddo ; enddo
  else
    ! Consistent element-quadrature assembly, matching the basal friction integration in CG_action
    ! (Lipscomb 2019 A4, eqs A26-A27): at each 2x2 Gauss point the nodal slope is interpolated with
    ! the bilinear basis, multiplied by the cell-mean thickness, and distributed back to the four
    ! cell-corner nodes weighted by that basis. Matching the friction's quadrature co-locates the
    ! driving-stress and friction grounding lines. The per-quadrature-point weight 1/4 |J_q| (CS%Jac)
    ! carries the element metric, exactly as CG_action's jac_wt = |J_q|/areaT does for friction; it
    ! reduces to 1/4 areaT on rectangular cells where |J_q| == areaT.
    ! Each element writes its four corner contributions into its own slot of taudx_b/taudy_b; the
    ! quadrature-point and four-element sums are then accumulated in diagonal pairs (exactly as
    ! CG_action assembles uret_b) so the load vector is bitwise invariant under 90 deg rotation.
    taudx_b(:,:,:) = 0.0 ; taudy_b(:,:,:) = 0.0
    do j=jsc-1,jec+1 ; do i=isc-1,iec+1
      if (ice_cell(i,j)) then
        He = max(ISS%h_shelf(i,j), CS%min_h_shelf)
        do jq=1,2 ; do iq=1,2
          qp = 2*(jq-1)+iq
          ! Per-quadrature-point area weight |J_q|/4, the same element-map metric CG_action applies
          ! to the friction integral (via jac_wt). On non-rectangular (e.g. lat/lon) cells this keeps
          ! the driving-stress and friction integrals on an identical metric; = 0.25*areaT when |J_q|==areaT.
          wq = 0.25 * CS%Jac(qp,i,j)
          ! Bilinear basis of the 4 corners at this quadrature point (same convention as CG_action).
          do jphi=1,2 ; do iphi=1,2
            ilq = 1 ; if (iq == iphi) ilq = 2
            jlq = 1 ; if (jq == jphi) jlq = 2
            phim(iphi,jphi) = xquad(ilq) * xquad(jlq)
          enddo ; enddo
          ! Interpolate the nodal slope to the quadrature point, summing corners in diagonal pairs.
          dsdx_qp = ((sx_n(i-1,j-1)*phim(1,1) + sx_n(i,j)*phim(2,2)) + &
                     (sx_n(i,j-1)*phim(2,1) + sx_n(i-1,j)*phim(1,2)))
          dsdy_qp = ((sy_n(i-1,j-1)*phim(1,1) + sy_n(i,j)*phim(2,2)) + &
                     (sy_n(i,j-1)*phim(2,1) + sy_n(i-1,j)*phim(1,2)))
          ! Per-corner contribution -rho g H grad(s) phi_m |J| w_q at this quadrature point.
          do jphi=1,2 ; do iphi=1,2
            txqp(iphi,jphi,qp) = -(((rho*grav)*He) * dsdx_qp) * (phim(iphi,jphi)*wq)
            tyqp(iphi,jphi,qp) = -(((rho*grav)*He) * dsdy_qp) * (phim(iphi,jphi)*wq)
          enddo ; enddo
        enddo ; enddo
        ! Sum the four quadrature points in diagonal pairs (qp 1+4, 2+3) into each corner's element slot.
        taudx_b(I-1,J-1,4) = (txqp(1,1,1)+txqp(1,1,4)) + (txqp(1,1,2)+txqp(1,1,3))
        taudx_b(I-1,J  ,2) = (txqp(1,2,1)+txqp(1,2,4)) + (txqp(1,2,2)+txqp(1,2,3))
        taudx_b(I  ,J-1,3) = (txqp(2,1,1)+txqp(2,1,4)) + (txqp(2,1,2)+txqp(2,1,3))
        taudx_b(I  ,J  ,1) = (txqp(2,2,1)+txqp(2,2,4)) + (txqp(2,2,2)+txqp(2,2,3))
        taudy_b(I-1,J-1,4) = (tyqp(1,1,1)+tyqp(1,1,4)) + (tyqp(1,1,2)+tyqp(1,1,3))
        taudy_b(I-1,J  ,2) = (tyqp(1,2,1)+tyqp(1,2,4)) + (tyqp(1,2,2)+tyqp(1,2,3))
        taudy_b(I  ,J-1,3) = (tyqp(2,1,1)+tyqp(2,1,4)) + (tyqp(2,1,2)+tyqp(2,1,3))
        taudy_b(I  ,J  ,1) = (tyqp(2,2,1)+tyqp(2,2,4)) + (tyqp(2,2,2)+tyqp(2,2,3))
      endif
    enddo ; enddo
    ! Assemble the four surrounding-element slots at each node in diagonal pairs (slots 1+4, 2+3).
    do J=jsc-1,jec ; do I=isc-1,iec
      taudx(I,J) = taudx(I,J) + ((taudx_b(I,J,1)+taudx_b(I,J,4)) + (taudx_b(I,J,2)+taudx_b(I,J,3)))
      taudy(I,J) = taudy(I,J) + ((taudy_b(I,J,1)+taudy_b(I,J,4)) + (taudy_b(I,J,2)+taudy_b(I,J,3)))
    enddo ; enddo
  endif

  ! Lateral-pressure (Neumann) boundary conditions at calving fronts and stress faces. This block
  ! is identical to calc_shelf_driving_stress; the front forcing is unchanged by the gradient scheme.
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then
      if (CS%ground_frac(i,j) == 1) then
        neumann_val = ((.5 * grav) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2 - &
                                      rhow * max(0.0, CS%bed_elev(i,j))**2))
      else
        neumann_val = (.5 * grav) * ((1-rho/rhow) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2))
      endif
      if ((CS%u_face_mask_bdry(I-1,j) == 2) .OR. &
        ((ISS%hmask(i-1,j) == 0 .OR. ISS%hmask(i-1,j) == 2) .AND. (CS%reentrant_x .OR. (i+i_off /= gisc)))) then
        taudx(I-1,J-1) = taudx(I-1,J-1) - .5 * G%dyCu(I-1,j) * neumann_val
        taudx(I-1,J) = taudx(I-1,J) - .5 * G%dyCu(I-1,j) * neumann_val
      endif
      if ((CS%u_face_mask_bdry(I,j) == 2) .OR. &
        ((ISS%hmask(i+1,j) == 0 .OR. ISS%hmask(i+1,j) == 2) .and. (CS%reentrant_x .OR. (i+i_off /= giec)))) then
        taudx(I,J-1) = taudx(I,J-1) + .5 * G%dyCu(I,j) * neumann_val
        taudx(I,J) = taudx(I,J) + .5 * G%dyCu(I,j) * neumann_val
      endif
      if ((CS%v_face_mask_bdry(i,J-1) == 2) .OR. &
        ((ISS%hmask(i,j-1) == 0 .OR. ISS%hmask(i,j-1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjsc)))) then
        taudy(I-1,J-1) = taudy(I-1,J-1) - .5 * G%dxCv(i,J-1) * neumann_val
        taudy(I,J-1) = taudy(I,J-1) - .5 * G%dxCv(i,J-1) * neumann_val
      endif
      if ((CS%v_face_mask_bdry(i,J) == 2) .OR. &
        ((ISS%hmask(i,j+1) == 0 .OR. ISS%hmask(i,j+1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjec)))) then
        taudy(I-1,J) = taudy(I-1,J) + .5 * G%dxCv(i,J) * neumann_val
        taudy(I,J) = taudy(I,J) + .5 * G%dxCv(i,J) * neumann_val
      endif
    endif
  enddo ; enddo

  ! Surface-slope diagnostic at cell centers: average the four corner-node slopes.
  if (CS%id_sx_shelf > 0 .or. CS%id_sy_shelf > 0 .or. CS%id_surf_slope_mag_shelf > 0) then
    do j=jsc,jec ; do i=isc,iec
      if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then
        CS%sx_shelf(i,j) = 0.25*((sx_n(I-1,J-1) + sx_n(I,J)) + (sx_n(I,J-1) + sx_n(I-1,J)))
        CS%sy_shelf(i,j) = 0.25*((sy_n(I-1,J-1) + sy_n(I,J)) + (sy_n(I,J-1) + sy_n(I-1,J)))
      else
        CS%sx_shelf(i,j) = 0.0 ; CS%sy_shelf(i,j) = 0.0
      endif
    enddo ; enddo
  endif

contains

  !> Option-3 (Lipscomb 2019) test for whether the edge between cells A=(ia,ja) and B=(ib,jb)
  !! contributes a surface-slope estimate: yes if both are ice-covered, or if an ice-covered cell
  !! lies higher in surface elevation than an ice-free land neighbor; no across ice/ice-free-ocean
  !! margins or where ice lies below ice-free land (a nunatak). An edge is also rejected if either
  !! cell lies outside a non-reentrant computational boundary, so the nodal gradient never differences
  !! across a solid wall into the (unphysical) across-wall halo. Without this guard the N/S (and E/W)
  !! wall nodes pick up a spurious cross-wall slope, breaking the meridional symmetry of channel
  !! configurations like MISMIP3D, matching the boundary handling in calc_shelf_driving_stress.
  !> True if both cells of an edge lie inside the global computational domain. An edge that fails
  !! this test has no CISM counterpart (CISM's halo BC is periodic, so its across-wall cells hold
  !! ice); it is excluded from both the numerator and the divisor, which is equivalent to a mirror
  !! boundary condition on the surface slope.
  logical function edge_in_domain(ia, ja, ib, jb)
    integer, intent(in) :: ia, ja, ib, jb
    edge_in_domain = .false.
    if (.not. CS%reentrant_x) then
      if ((ia+i_off < gisc) .or. (ia+i_off > giec) .or. &
          (ib+i_off < gisc) .or. (ib+i_off > giec)) return
    endif
    if (.not. CS%reentrant_y) then
      if ((ja+j_off < gjsc) .or. (ja+j_off > gjec) .or. &
          (jb+j_off < gjsc) .or. (jb+j_off > gjec)) return
    endif
    edge_in_domain = .true.
  end function edge_in_domain

  !> True if an in-domain edge carries a meaningful surface gradient: ice in both cells, or ice
  !! standing above an ice-free land cell (CISM HO_GRADIENT_MARGIN_HYBRID). An edge that fails this
  !! test is the option-3 margin case: it contributes nothing but still counts in the divisor, so
  !! the surviving parallel edge is halved, as in CISM's centered edge-gradient average.
  logical function edge_ice_ok(ia, ja, ib, jb)
    integer, intent(in) :: ia, ja, ib, jb
    edge_ice_ok = .false.
    if (ice_cell(ia,ja) .and. ice_cell(ib,jb)) then
      edge_ice_ok = .true.
    elseif (ice_cell(ia,ja) .and. land_cell(ib,jb)) then
      edge_ice_ok = (S(ia,ja) > S(ib,jb))
    elseif (ice_cell(ib,jb) .and. land_cell(ia,ja)) then
      edge_ice_ok = (S(ib,jb) > S(ia,ja))
    endif
  end function edge_ice_ok

  !> Metric-exact lumped (row-sum) nodal mass contribution of cell (ic,jc) to the node at its
  !! (icorner,jcorner) corner: h * integral of that corner's bilinear basis against the element
  !! Jacobian, (h/4) * sum_qp phi_corner(qp) * Jac(qp). Two-point Gauss is exact for this integrand
  !! on the bilinear locally-orthogonal map (phi*Jac is at most quadratic per direction), so the
  !! result carries the lat/lon metric exactly, matching the per-quadrature-point |J_q| weighting
  !! CG_action uses. Reduces to 0.25*areaT*h on rectangular cells where Jac(qp) == areaT.
  real function lumped_corner_mass(ic, jc, icorner, jcorner)
    integer, intent(in) :: ic, jc        !< Cell indices of the contributing cell
    integer, intent(in) :: icorner, jcorner !< The cell corner (1=W/S, 2=E/N) that is the node
    real :: pj(4)  ! phi_corner(qp) * Jac(qp) at the four quadrature points [L2 ~> m2]
    integer :: iq2, jq2, qq, il, jl
    do jq2=1,2 ; do iq2=1,2
      qq = 2*(jq2-1)+iq2
      il = 1 ; if (iq2 == icorner) il = 2
      jl = 1 ; if (jq2 == jcorner) jl = 2
      pj(qq) = (xquad(il)*xquad(jl)) * CS%Jac(qq,ic,jc)
    enddo ; enddo
    ! Gauss weight 1/4 per QP, summed in diagonal pairs (1+4)+(2+3) to match CG_action's
    ! rotation-invariant quadrature; times the cell-mean thickness.
    lumped_corner_mass = (0.25 * ((pj(1)+pj(4)) + (pj(2)+pj(3)))) * &
                         max(ISS%h_shelf(ic,jc), CS%min_h_shelf)
  end function lumped_corner_mass

end subroutine calc_shelf_driving_stress_vertex

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
  real :: gl_w_qp        ! Quadrant-GLP grounded weight at the quadrature point, bilinear-interpolated
                         ! from CS%f_ground_node; used only when CS%gl_quad_friction [nondim]
  real :: drag_newt_qp   ! Newton basal drag coefficient at quadrature point [R Z T-1 ~> kg m-2 s-1]
  real :: inner_dot_qp   ! u^k_qp · δu_qp inner product for Newton basal drag [L2 T-2 ~> m2 s-2]
  real :: bcoef_loc, dnewt_loc ! Local (nodal-diagonal) basal Picard drag [R L2 Z T-1 ~> kg s-1] and
                         ! Newton tangent factor [R Z T ~> kg m-2 s] at a node (CISM_FRICTION='local')
  real :: idot_loc       ! u^k_node · δu_node inner product for the local Newton drag [L2 T-2 ~> m2 s-2]
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
  logical :: fv_sub_fric ! Local flag for I2N_FRICTION (corner H and fls fields)
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
  fv_sub_fric = CS%i2n_friction
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
        if (CS%gl_quad_friction) then
          ! Quadrant GLP: the element friction path handles every cell with any grounded area,
          ! and the per-quadrature-point weight gl_w_qp (below) carries the sub-cell structure.
          grounded_qp = CS%f_ground_cell(i,j) > 0.0
        else
          grounded_qp = merge(CS%ground_frac(i,j) >= 1.0, CS%ground_frac(i,j) > 0.0, CS%GL_regularize)
        endif
        if (grounded_qp) then
          ! DG: h_gp is only used for flotation and the effective pressure.
          if (do_DG) then
            h_gp = ((CS%h_nodal(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                    (CS%h_nodal(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
                   ((CS%h_nodal(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                    (CS%h_nodal(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
            h_gp = max(h_gp, CS%min_h_shelf)
            bed_gp = ((CS%bed_node(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                      (CS%bed_node(I,J)     * (xquad(iq)   * xquad(jq))))  + &
                     ((CS%bed_node(I,J-1)   * (xquad(iq)   * xquad(3-jq))) + &
                      (CS%bed_node(I-1,J)   * (xquad(3-iq) * xquad(jq))))
            ! Under quadrant GLP the smooth nodal weight gl_w_qp sets the grounded contribution,
            ! so the binary per-quadrature-point test must not gate the friction; it still selects
            ! a physical (grounded-only) effective pressure for the Coulomb law.
            if (.not. CS%gl_quad_friction) grounded_qp = (dens_ratio * h_gp - bed_gp > 0)
            if ((dens_ratio * h_gp - bed_gp > 0) .and. CS%CoulombFriction) then
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
          ! Quadrant-GLP grounded weight at this quadrature point: bilinear interpolation of the
          ! nodal grounded fraction (same corner basis as h_gp above). Harmless when the toggle is
          ! off (f_ground_node is zero and the merge below discards it).
          gl_w_qp = ((CS%f_ground_node(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                     (CS%f_ground_node(I,J)     * (xquad(iq)   * xquad(jq))))  + &
                    ((CS%f_ground_node(I,J-1)   * (xquad(iq)   * xquad(3-jq))) + &
                     (CS%f_ground_node(I-1,J)   * (xquad(3-iq) * xquad(jq))))
          ! Apply ground fraction scaling (replaces external scaling of basal_traction).
          ! Under GL_regularize, GL cells (0 < ground_frac < 1) get sub-grid-aware basal
          ! handling via the GL branch below, so the cell-level scaling must collapse to 1.0
          ! there to avoid double-counting (matches pre-refactor behavior when ground_frac
          ! was forced to 1.0 in GL cells). Under quadrant GLP the nodal weight gl_w_qp is used
          ! at every grounded cell instead.
          basal_coef_qp = basal_coef_qp * &
              merge(gl_w_qp, merge(1.0, CS%ground_frac(i,j), CS%GL_regularize), CS%gl_quad_friction)
          if (use_newton) then
            drag_newt_qp = drag_newt_qp * &
                merge(gl_w_qp, merge(1.0, CS%ground_frac(i,j), CS%GL_regularize), CS%gl_quad_friction)
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

          if (grounded_qp .and. .not. CS%local_basal_friction) then
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

      if (CS%GL_regularize .and. .not. CS%gl_quad_friction .and. .not. CS%local_basal_friction .and. &
          CS%ground_frac(i,j) > 0.0 .and. CS%ground_frac(i,j) < 1.0) then
        ! Subgrid grounding-line: evaluate basal friction at each grounded sub-quadrature point.
        ! Picard and Newton Jacobian are both computed inside CG_action_subgrid_basal.
        Hcell(:,:) = H_node(I-1:I,J-1:J)
        if (CS%use_sep2) then
          if (fv_sub_fric) then
            call CG_action_sep2_basal(CS, G, US, Hcell, &
                u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
                u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
                bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
                G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
                h_nodal_cell=CS%H_corner(I-1:I,J-1:J), &
                fls_cell=CS%fls_corner(I-1:I,J-1:J))
          elseif (do_DG) then
            call CG_action_sep2_basal(CS, G, US, Hcell, &
                u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
                u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
                bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
                G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
                use_DG=.true., h_nodal_cell=CS%h_nodal(i,j,:,:), &
                bed_corners=CS%bed_node(I-1:I,J-1:J))
          else
            call CG_action_sep2_basal(CS, G, US, Hcell, &
                u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
                u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
                bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
                G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j))
          endif
        elseif (fv_sub_fric) then
          call CG_action_subgrid_basal(CS, G, US, Phisub, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
              bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
              h_nodal_cell=CS%H_corner(I-1:I,J-1:J), &
              fls_cell=CS%fls_corner(I-1:I,J-1:J))
        elseif (do_DG) then
          ! h_nodal_cell is only used for flotation and the effective pressure.
          call CG_action_subgrid_basal(CS, G, US, Phisub, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              u_shlf(I-1:I,J-1:J), v_shlf(I-1:I,J-1:J), &
              bathyT(i,j), dens_ratio, i, j, fB_e, use_newton, Usub, Vsub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
              use_DG=.true., &
              h_nodal_cell=CS%h_nodal(i,j,:,:), &
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

  ! Local (nodal-diagonal) basal friction (CISM HO_ASSEMBLE_BETA_LOCAL): each node's drag uses only
  ! its own velocity, beta, and grounded fraction f_ground_node -- no element integration or neighbor
  ! coupling. Replaces the per-quadrature-point friction gated off above.
  if (CS%local_basal_friction) then
    do J=js-1,je ; do I=is-1,ie
      if (CS%f_ground_node(I,J) > 0.0) then
        call compute_basal_coef_node(CS, G, US, I, J, u_curr(I,J), v_curr(I,J), use_newton, &
                                     bcoef_loc, dnewt_loc)
        idot_loc = (u_curr(I,J)*u_shlf(I,J)) + (v_curr(I,J)*v_shlf(I,J))
        if (umask(I,J) == 1) uret(I,J) = uret(I,J) + &
            ((bcoef_loc * u_shlf(I,J)) + (dnewt_loc * u_curr(I,J) * idot_loc))
        if (vmask(I,J) == 1) vret(I,J) = vret(I,J) + &
            ((bcoef_loc * v_shlf(I,J)) + (dnewt_loc * v_curr(I,J) * idot_loc))
      endif
    enddo ; enddo
  endif

end subroutine CG_action

!> Compute subgrid grounding-line basal traction nodal contributions for a CG action.
!! Evaluates basal friction (Picard and Newton Jacobian) at each grounded sub-quadrature point.
!! The sub-qp flotation test accounts for partial grounding; no external ground_frac scaling needed.
subroutine CG_action_subgrid_basal(CS, G, US, Phisub, H, U_curr, V_curr, U_delta, V_delta, &
                                   bathyT, dens_ratio, i_elem, j_elem, fB_e, use_newton, Ucontr, Vcontr, &
                                   dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT, &
                                   use_DG, h_nodal_cell, bed_corners, fls_cell)
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
  real, dimension(2,2), optional, intent(in) :: h_nodal_cell !< Q1 thickness at the 4 cell corners, used
                                          !! only as a flotation measure (sub-qp gate + effective
                                          !! pressure) [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: bed_corners !< Bed elevation at element corners [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: fls_cell !< Flotation deficit r*h - bed at the 4 cell
                                              !! corners (I2N_FRICTION). When present the
                                              !! sub-point flotation test and the effective pressure
                                              !! are taken from this field and h_nodal_cell supplies
                                              !! the matching corner thickness [Z ~> m]

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
  logical :: active_qp   ! Whether this sub-qp contributes basal traction
  logical :: do_fvsub    ! Local flag for the FV sub-element mode (fls_cell supplied)
  real :: fls_loc        ! Flotation deficit r*h - bed at sub-qp [Z ~> m]
  real :: rho_ocean_g_LtoZ ! US%L_to_Z * density_ocean_avg * g_Earth [R L Z-1 T-2]
  integer :: nsub, i, j, qx, qy, m, n

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  coef_prefactor = CS%coef_prefactor(i_elem,j_elem)
  min_trac_area  = CS%min_basal_traction * G%areaT(i_elem,j_elem)
  eps_vel2 = CS%eps_glen_min**2 * ((G%dxT(i_elem,j_elem)**2) + (G%dyT(i_elem,j_elem)**2))


  do_DG = .false.
  if (present(use_DG)) do_DG = use_DG
  do_fvsub = present(fls_cell)
  if (do_fvsub) then
    rho_ocean_g_LtoZ = US%L_to_Z * CS%density_ocean_avg * CS%g_Earth
  elseif (do_DG) then
    rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
    rho_ice_g_LtoZ = US%L_to_Z * CS%density_ice * CS%g_Earth
  endif

  Ucontr_sub(:,:,:,:) = 0.0 ; Vcontr_sub(:,:,:,:) = 0.0

  do j=1,nsub ; do i=1,nsub
    U_qp_nd(:,:,:,:) = 0.0 ; V_qp_nd(:,:,:,:) = 0.0
    do qy=1,2 ; do qx=1,2
      if (do_fvsub) then
        ! FV sub-element mode: thickness and flotation deficit share one interpolant (the same
        ! bilinear Phisub weights that the SEP3 sub-point flotation test uses), so the grounding
        ! line, the effective pressure and the surface kink all key off one field.
        hloc = ((Phisub(qx,qy,i,j,1,1)*h_nodal_cell(1,1)) + (Phisub(qx,qy,i,j,2,2)*h_nodal_cell(2,2))) + &
               ((Phisub(qx,qy,i,j,1,2)*h_nodal_cell(1,2)) + (Phisub(qx,qy,i,j,2,1)*h_nodal_cell(2,1)))
        fls_loc = ((Phisub(qx,qy,i,j,1,1)*fls_cell(1,1)) + (Phisub(qx,qy,i,j,2,2)*fls_cell(2,2))) + &
                  ((Phisub(qx,qy,i,j,1,2)*fls_cell(1,2)) + (Phisub(qx,qy,i,j,2,1)*fls_cell(2,1)))
        bed_sub = bathyT  ! unused; the bed is implicit in fls
      elseif (do_DG) then
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

      ! Plain flotation test at the sub-quadrature point.
      if (do_fvsub) then
        active_qp = (fls_loc > 0)
      else
        active_qp = (dens_ratio * hloc - bed_sub > 0)
      endif
      if (active_qp) then  ! grounded sub-qp
        u_curr_loc  = (((Phisub(qx,qy,i,j,1,1)*U_curr(1,1))  + (Phisub(qx,qy,i,j,2,2)*U_curr(2,2)))  + &
                       ((Phisub(qx,qy,i,j,1,2)*U_curr(1,2))  + (Phisub(qx,qy,i,j,2,1)*U_curr(2,1))))
        v_curr_loc  = (((Phisub(qx,qy,i,j,1,1)*V_curr(1,1))  + (Phisub(qx,qy,i,j,2,2)*V_curr(2,2)))  + &
                       ((Phisub(qx,qy,i,j,1,2)*V_curr(1,2))  + (Phisub(qx,qy,i,j,2,1)*V_curr(2,1))))
        u_delta_loc = (((Phisub(qx,qy,i,j,1,1)*U_delta(1,1)) + (Phisub(qx,qy,i,j,2,2)*U_delta(2,2))) + &
                       ((Phisub(qx,qy,i,j,1,2)*U_delta(1,2)) + (Phisub(qx,qy,i,j,2,1)*U_delta(2,1))))
        v_delta_loc = (((Phisub(qx,qy,i,j,1,1)*V_delta(1,1)) + (Phisub(qx,qy,i,j,2,2)*V_delta(2,2))) + &
                       ((Phisub(qx,qy,i,j,1,2)*V_delta(1,2)) + (Phisub(qx,qy,i,j,2,1)*V_delta(2,1))))

        unorm2_loc = ((u_curr_loc**2) + (v_curr_loc**2)) + eps_vel2

        ! Compute Coulomb fB at this sub-qp when the effective pressure varies within the cell
        if (do_fvsub .and. CS%CoulombFriction) then
          fB_local = compute_fB_from_N( &
              subgrid_effective_pressure(fls_loc, hloc, dens_ratio, rho_ocean_g_LtoZ), &
              CS%C_basal_friction(i_elem,j_elem), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
              CS%CF_PostPeak, CS%n_basal_fric)
        elseif (do_DG .and. CS%CoulombFriction) then
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
  elseif (CoulombFriction .and. (fB_e < 0.0)) then
    ! Zero effective pressure (FB_NO_COULOMB_DRAG): the Coulomb law gives exactly zero basal drag
    ! and its Newton tangent vanishes with it, since the drag is identically zero in u. Taking this
    ! branch avoids forming the divergent fB expression at all. min_trac_area reproduces what the
    ! Coulomb branch below would do with raw_coef = 0: the floor when one is set, zero otherwise.
    basal_coef = min_trac_area  ;  drag_newt = 0.0
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

!> Local (nodal) basal drag coefficient and Newton tangent factor at B-grid node (I,J) for
!! CISM_FRICTION='local' (CISM HO_ASSEMBLE_BETA_LOCAL). beta is built from the pre-computed nodal
!! prefactor (areaBu*C_node) and nodal Coulomb fB, the nodal velocity magnitude, then scaled by the
!! nodal grounded fraction f_ground_node. The drag is purely diagonal: tau_b at the node = bcoef*u.
subroutine compute_basal_coef_node(CS, G, US, I, J, u_c, v_c, use_newton, bcoef, dnewt)
  type(ice_shelf_dyn_CS), intent(in) :: CS  !< Ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G   !< The grid structure
  type(unit_scale_type),  intent(in) :: US  !< Unit conversion factors
  integer, intent(in) :: I, J               !< B-grid node indices
  real,    intent(in) :: u_c, v_c           !< Current nodal velocity components [L T-1 ~> m s-1]
  logical, intent(in) :: use_newton         !< If true, evaluate the Newton tangent factor
  real,    intent(out) :: bcoef             !< Picard diagonal drag at the node [R L2 Z T-1 ~> kg s-1]
  real,    intent(out) :: dnewt             !< Newton drag tangent factor [R Z T ~> kg m-2 s]; 0 without Newton

  real :: eps2     ! Velocity regularization squared at the node [L2 T-2 ~> m2 s-2]
  real :: mintrac  ! min_basal_traction * area_node floor at the node [R L2 Z T-1 ~> kg s-1]
  real :: unorm2   ! Regularized |u|^2 at the node [L2 T-2 ~> m2 s-2]

  ! Scale the regularization and traction floor by the same nodal control volume (CS%area_node) the
  ! drag prefactor uses, not areaBu, so every term in the local node balance shares one control volume
  ! and the domain-edge nodes stay consistent with the interior (meridional symmetry).
  eps2    = CS%eps_glen_min**2 * CS%area_node(I,J)
  mintrac = CS%min_basal_traction * CS%area_node(I,J)
  unorm2  = ((u_c**2) + (v_c**2)) + eps2
  ! CISM HO_BETA_LIMIT_ABSOLUTE (its default): scale by the grounded fraction first, then floor, so
  ! that a node with any grounded area keeps at least MIN_BASAL_TRACTION. Where the floor binds the
  ! drag is constant in |u|, so the Newton tangent factor vanishes, as in compute_basal_coef.
  call compute_basal_coef(unorm2, CS%coef_prefactor_node(I,J), 0.0, CS%fB_node(I,J), &
      CS%n_basal_fric, CS%CoulombFriction, CS%CF_PostPeak, US%L_T_to_m_s, use_newton, bcoef, dnewt)
  bcoef = bcoef * CS%f_ground_node(I,J)
  dnewt = dnewt * CS%f_ground_node(I,J)
  if ((CS%f_ground_node(I,J) > 0.0) .and. (bcoef < mintrac)) then
    bcoef = mintrac ; dnewt = 0.0
  endif
end subroutine compute_basal_coef_node

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

  real :: fN  ! Effective pressure [R Z L T-2 ~> Pa]

  ! N is already floored at CF_MinN by coulomb_effective_pressure, so no further floor is applied
  ! here; compute_fB_from_N returns FB_NO_COULOMB_DRAG if that leaves it at zero.
  fN = coulomb_effective_pressure(h_local, bed_local, rho_oi_ratio, rho_ice_g_LtoZ, CF_MinN)
  compute_fB_local = compute_fB_from_N(fN, C_basal, alpha_coulomb, CF_Max, 0.0, CF_PostPeak, n_basal_fric)
end function compute_fB_local

!> The Coulomb effective pressure N = rho_i*g*(H - H_f) at a point, with H_f the flotation
!! thickness, floored at N_min and capped at the overburden pressure rho_i*g*H. This is the
!! p_ocean_penetration = 1 case of the Leguy et al. (2021) effective pressure, and matches CISM's
!! calc_effective_pressure (glissade_basal_traction) with that exponent, including its cap of N to
!! [0, overburden]. The upper cap is redundant here (H_f >= 0 already gives N <= rho_i*g*H) but is
!! kept explicit so the correspondence with CISM is visible. It is applied before the floor, so that
!! passing N_min = CF_MinN reproduces the unclamped expression exactly for any thickness.
pure real function coulomb_effective_pressure(h_local, bed_local, rho_oi_ratio, rho_ice_g_LtoZ, N_min)
  real, intent(in) :: h_local        !< Ice thickness at the point [Z ~> m]
  real, intent(in) :: bed_local      !< Bed elevation (positive below sea level) at the point [Z ~> m]
  real, intent(in) :: rho_oi_ratio   !< density_ocean_avg / density_ice [nondim]
  real, intent(in) :: rho_ice_g_LtoZ !< US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  real, intent(in) :: N_min          !< Minimum effective pressure [R Z L T-2 ~> Pa]

  real :: Hf  ! Flotation thickness [Z ~> m]

  Hf = max(rho_oi_ratio * bed_local, 0.0)
  coulomb_effective_pressure = max(min(rho_ice_g_LtoZ * (h_local - Hf), rho_ice_g_LtoZ * h_local), N_min)
end function coulomb_effective_pressure


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
  real :: gl_w_qp        ! Quadrant-GLP grounded weight at the quadrature point, bilinear-interpolated
                         ! from CS%f_ground_node; used only when CS%gl_quad_friction [nondim]
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
  logical :: fv_sub_fric ! Local flag for I2N_FRICTION (corner H and fls fields)
  logical :: grounded_qp ! Whether this quadrature point is grounded
  real :: bcoef_loc, dnewt_loc ! Local (nodal-diagonal) basal Picard drag [R L2 Z T-1 ~> kg s-1] and
                         ! Newton tangent factor [R Z T ~> kg m-2 s] at a node (CISM_FRICTION='local')
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
  fv_sub_fric = CS%i2n_friction
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
      if (CS%gl_quad_friction) then
        ! Mirror the CG_action operator: element friction path for any grounded cell, weighted
        ! per quadrature point by gl_w_qp so the preconditioner diagonal matches the operator.
        grounded_qp = CS%f_ground_cell(i,j) > 0.0
      else
        grounded_qp = merge(CS%ground_frac(i,j) >= 1.0, CS%ground_frac(i,j) > 0.0, CS%GL_regularize)
      endif
      if (grounded_qp) then
        if (do_DG) then
          ! h_gp is used only as a flotation measure (gate + effective pressure), so it
          ! is read from the cell's own nodal thickness.
          h_gp = ((CS%h_nodal(i,j,1,1) * (xquad(3-iq) * xquad(3-jq))) + &
                  (CS%h_nodal(i,j,2,2) * (xquad(iq)   * xquad(jq))))  + &
                 ((CS%h_nodal(i,j,2,1) * (xquad(iq)   * xquad(3-jq))) + &
                  (CS%h_nodal(i,j,1,2) * (xquad(3-iq) * xquad(jq))))
          h_gp = max(h_gp, CS%min_h_shelf)
          ! Bed at the QP with the same rotation-paired bilinear pattern as
          ! u_curr_qp below, for symmetric rotation-cancel structure.
          bed_gp = ((CS%bed_node(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                    (CS%bed_node(I,J)     * (xquad(iq)   * xquad(jq))))  + &
                   ((CS%bed_node(I,J-1)   * (xquad(iq)   * xquad(3-jq))) + &
                    (CS%bed_node(I-1,J)   * (xquad(3-iq) * xquad(jq))))
          ! See CG_action: under quadrant GLP the smooth nodal weight sets the grounded
          ! contribution, so the binary per-QP test must not gate the diagonal; it still selects
          ! a physical (grounded-only) effective pressure for the Coulomb law.
          if (.not. CS%gl_quad_friction) grounded_qp = (dens_ratio * h_gp - bed_gp > 0)
          if ((dens_ratio * h_gp - bed_gp > 0) .and. CS%CoulombFriction) then
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
        ! Quadrant-GLP grounded weight at this QP (same corner basis as h_gp); matches CG_action.
        gl_w_qp = ((CS%f_ground_node(I-1,J-1) * (xquad(3-iq) * xquad(3-jq))) + &
                   (CS%f_ground_node(I,J)     * (xquad(iq)   * xquad(jq))))  + &
                  ((CS%f_ground_node(I,J-1)   * (xquad(iq)   * xquad(3-jq))) + &
                   (CS%f_ground_node(I-1,J)   * (xquad(3-iq) * xquad(jq))))
        basal_coef_qp = basal_coef_qp * &
            merge(gl_w_qp, merge(1.0, CS%ground_frac(i,j), CS%GL_regularize), CS%gl_quad_friction)
        drag_newt_qp  = drag_newt_qp  * &
            merge(gl_w_qp, merge(1.0, CS%ground_frac(i,j), CS%GL_regularize), CS%gl_quad_friction)
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

          if (grounded_qp .and. .not. CS%local_basal_friction) then
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

          if (grounded_qp .and. .not. CS%local_basal_friction) then
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

    if (CS%GL_regularize .and. .not. CS%gl_quad_friction .and. .not. CS%local_basal_friction .and. &
        CS%ground_frac(i,j) > 0.0 .and. CS%ground_frac(i,j) < 1.0) then
      ! Subgrid grounding-line: evaluate basal friction diagonal at each grounded sub-quadrature point.
      ! Returns separate u_diag_sub and v_diag_sub (differ in Newton term: u^2 vs v^2).
      ! The sub-qp flotation test handles grounding fraction; no external ground_frac scaling needed.
      Hcell(:,:) = H_node(I-1:I,J-1:J)
      if (CS%use_sep2) then
        if (fv_sub_fric) then
          call CG_diagonal_sep2_basal(CS, G, US, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
              h_nodal_cell=CS%H_corner(I-1:I,J-1:J), &
              fls_cell=CS%fls_corner(I-1:I,J-1:J))
        elseif (do_DG) then
          call CG_diagonal_sep2_basal(CS, G, US, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
              use_DG=.true., h_nodal_cell=CS%h_nodal(i,j,:,:), &
              bed_corners=CS%bed_node(I-1:I,J-1:J))
        else
          call CG_diagonal_sep2_basal(CS, G, US, Hcell, &
              u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
              CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
              G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j))
        endif
      elseif (fv_sub_fric) then
        call CG_diagonal_subgrid_basal(CS, G, US, Phisub, Hcell, &
            u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
            CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
            G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
            h_nodal_cell=CS%H_corner(I-1:I,J-1:J), &
            fls_cell=CS%fls_corner(I-1:I,J-1:J))
      elseif (do_DG) then
        ! h_nodal_cell is only used for flotation and the effective pressure.
        call CG_diagonal_subgrid_basal(CS, G, US, Phisub, Hcell, &
            u_curr(I-1:I,J-1:J), v_curr(I-1:I,J-1:J), &
            CS%bed_elev(i,j), dens_ratio, i, j, fB_e, u_diag_sub, v_diag_sub, &
            G%dxCv(i,j-1), G%dxCv(i,j), G%dyCu(i-1,j), G%dyCu(i,j), G%IareaT(i,j), &
            use_DG=.true., &
            h_nodal_cell=CS%h_nodal(i,j,:,:), &
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

  ! Local (nodal-diagonal) basal friction (CISM HO_ASSEMBLE_BETA_LOCAL): the diagonal of the nodal
  ! drag bcoef*u (+ Newton tangent dnewt*u^2), matching the term added in CG_action.
  if (CS%local_basal_friction) then
    do J=jsc-2,jec+1 ; do I=isc-2,iec+1
      if (CS%f_ground_node(I,J) > 0.0) then
        call compute_basal_coef_node(CS, G, US, I, J, u_curr(I,J), v_curr(I,J), CS%doing_newton, &
                                     bcoef_loc, dnewt_loc)
        if (CS%umask(I,J) == 1) u_diagonal(I,J) = u_diagonal(I,J) + (bcoef_loc + (dnewt_loc * u_curr(I,J)**2))
        if (CS%vmask(I,J) == 1) v_diagonal(I,J) = v_diagonal(I,J) + (bcoef_loc + (dnewt_loc * v_curr(I,J)**2))
      endif
    enddo ; enddo
  endif

end subroutine matrix_diagonal

!> Compute subgrid grounding-line basal traction contributions for the preconditioner diagonal.
!! Evaluates friction at each grounded sub-quadrature point. Returns separate u and v diagonals
!! because the Newton term uses u^2 for the u-block and v^2 for the v-block.
!! The sub-qp flotation test handles partial grounding; no external ground_frac scaling needed.
subroutine CG_diagonal_subgrid_basal(CS, G, US, Phisub, H_node, U_curr, V_curr, &
                                     bathyT, dens_ratio, i_elem, j_elem, fB_e, u_diag, v_diag, &
                                     dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT, &
                                     use_DG, h_nodal_cell, bed_corners, fls_cell)
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
  real, dimension(2,2), optional, intent(in) :: h_nodal_cell !< Q1 thickness at the 4 cell corners, used
                                          !! only as a flotation measure (sub-qp gate + effective
                                          !! pressure) [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: bed_corners !< Bed elevation at element corners [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: fls_cell !< Flotation deficit r*h - bed at the 4 cell
                                              !! corners (I2N_FRICTION). When present the
                                              !! sub-point flotation test and the effective pressure
                                              !! are taken from this field and h_nodal_cell supplies
                                              !! the matching corner thickness [Z ~> m]

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
  logical :: active_qp   ! Whether this sub-qp contributes basal traction
  logical :: do_fvsub    ! Local flag for the FV sub-element mode (fls_cell supplied)
  real :: fls_loc        ! Flotation deficit r*h - bed at sub-qp [Z ~> m]
  real :: rho_ocean_g_LtoZ ! US%L_to_Z * density_ocean_avg * g_Earth [R L Z-1 T-2]
  integer :: nsub, i, j, qx, qy, m, n

  nsub    = size(Phisub, 3)
  subarea = 1.0 / real(nsub)**2

  coef_prefactor = CS%coef_prefactor(i_elem,j_elem)
  min_trac_area  = CS%min_basal_traction * G%areaT(i_elem,j_elem)
  eps_vel2 = CS%eps_glen_min**2 * ((G%dxT(i_elem,j_elem)**2) + (G%dyT(i_elem,j_elem)**2))


  do_DG = .false.
  if (present(use_DG)) do_DG = use_DG
  do_fvsub = present(fls_cell)
  if (do_fvsub) then
    rho_ocean_g_LtoZ = US%L_to_Z * CS%density_ocean_avg * CS%g_Earth
  elseif (do_DG) then
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

      ! Plain flotation test at the sub-quadrature point (matches CG_action_subgrid_basal
      ! so the preconditioner diagonal stays consistent with the residual).
      if (do_fvsub) then
        active_qp = (fls_loc > 0)
      else
        active_qp = (dens_ratio * hloc - bed_sub > 0)
      endif
      if (active_qp) then  ! grounded sub-qp
        u_curr_loc = (((Phisub(qx,qy,i,j,1,1)*U_curr(1,1)) + (Phisub(qx,qy,i,j,2,2)*U_curr(2,2))) + &
                      ((Phisub(qx,qy,i,j,1,2)*U_curr(1,2)) + (Phisub(qx,qy,i,j,2,1)*U_curr(2,1))))
        v_curr_loc = (((Phisub(qx,qy,i,j,1,1)*V_curr(1,1)) + (Phisub(qx,qy,i,j,2,2)*V_curr(2,2))) + &
                      ((Phisub(qx,qy,i,j,1,2)*V_curr(1,2)) + (Phisub(qx,qy,i,j,2,1)*V_curr(2,1))))

        unorm2_loc = ((u_curr_loc**2) + (v_curr_loc**2)) + eps_vel2

        ! Compute Coulomb fB at this sub-qp when the effective pressure varies within the cell
        if (do_fvsub .and. CS%CoulombFriction) then
          fB_local = compute_fB_from_N( &
              subgrid_effective_pressure(fls_loc, hloc, dens_ratio, rho_ocean_g_LtoZ), &
              CS%C_basal_friction(i_elem,j_elem), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
              CS%CF_PostPeak, CS%n_basal_fric)
        elseif (do_DG .and. CS%CoulombFriction) then
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

!> Build the SEP2 sub-element quadrature of one grounding-line cell (Seroussi et al. 2014
!! SEP2, extended to quadrilaterals). The cell is fanned into 4 triangles at its center;
!! the bilinear flotation deficit f is linear on each, so the grounding line is a straight
!! cut separating the unique minority-sign vertex (triangle piece) from the other two
!! (quad piece, collapsed exactly when the cut passes through a vertex). QPs carry parent
!! Q1 corner-basis weights (beta) and reference-space measures (wref); no physical metric
!! enters here. All formulas are keyed to vertex roles and grouped in symmetry orbits so
!! outputs are bitwise-covariant under grid rotations and reflections.
!!
!! KNOWN LIMITATION: assembled basis-weighted quantities are discontinuous in the sub-cell
!! grounding-line position, although the partition geometry is not. The two rules below are
!! each exact for the parent Jacobian, so the weight sums -- hence areas and ground_frac --
!! vary continuously through every topology change. They are not exact for beta_i*beta_j*J,
!! which is total degree 4 (degree 5 per variable on a quad piece after the bilinear sub-map),
!! while the 3-pt triangle rule is exact to degree 2 and the 2x2 tensor rule to degree 3 per
!! variable. So the same region integrated by different rules returns different values, and the
!! assembled basal friction and driving stress step whenever a branch change re-assigns the
!! rules. The largest case is the centre deficit f_C changing sign, which switches all four
!! parent triangles at once: two go uncut <-> centre-cut, and two swap the minority vertex
!! between the corner-cut branches, which exchanges the triangle and quad rules over the same
!! two pieces. Measured on MISMIP+ at 10 km (dg_thin/runs/glsweep.sh, z18_glsmooth.py): a 2.5%
!! step in the local velocity at one grounding-line position, against a 0.02% background, not
!! resolved by refining the sweep 5x, and absent under GROUNDING_LINE_SUBGRID_SCHEME="SEP3".
!!
!! Raising the quadrature to 3x3 on quad pieces (exact to degree 5 per variable) and a
!! degree-4-exact rule on triangle pieces would remove it: rules that are exact for the
!! integrand agree whatever the parameterization, so a branch change cannot move the answer.
!! Using one rule family everywhere is NOT sufficient on its own -- which corner of a collapsed
!! quad is doubled up follows the branch's (X,Y,Z) role assignment, not the piece's geometry,
!! so the same piece is parameterized differently in the two branches that produce it.
!! Not done because it is not the cause of the steady grounding-line wobble (SEP3 removes these
!! steps and makes the wobble 3-8x worse), P75R reversibility already passes as this stands, and
!! the measured sensitivity to the quadrature choice away from a transition is about 0.2%.
subroutine sep2_cell_qps(f, nqp, beta, wref, qp_grounded)
  real, dimension(4),     intent(in)  :: f    !< Flotation deficit r*h - bed at the cell corners,
                                              !! ordered SW, SE, NW, NE [Z ~> m]
  integer, dimension(4),  intent(out) :: nqp  !< Number of QPs per parent triangle (3 or 7)
  real, dimension(4,7,4), intent(out) :: beta !< Corner-basis weights per (corner, QP, triangle),
                                              !! triangles ordered S, E, N, W [nondim]
  real, dimension(7,4),   intent(out) :: wref !< Reference-space measure per (QP, triangle) [nondim]
  logical, dimension(7,4), intent(out) :: qp_grounded !< True where the QP lies in a grounded piece

  ! Base-corner ids (A,B) of triangles S,E,N,W in counterclockwise order.
  integer, dimension(4), parameter :: iA = (/ 1, 2, 4, 3 /) ! First (CCW) base corner per triangle
  integer, dimension(4), parameter :: iB = (/ 2, 4, 3, 1 /) ! Second (CCW) base corner per triangle
  real, dimension(4) :: bC, bA, bB ! Corner-basis weights of the role vertices [nondim]
  real :: fC              ! Deficit at the cell center (bilinear value = corner mean) [Z ~> m]
  logical :: gC, gA, gB   ! Tie-broken grounded states (f > 0) of the role vertices
  integer :: k

  ! Diagonal-pair grouping: the two diagonals map to each other under any rotation/reflection.
  fC = 0.5 * ((0.5 * (f(1) + f(4))) + (0.5 * (f(2) + f(3))))
  gC = (fC > 0.0)
  bC(:) = 0.25

  do k=1,4
    gA = (f(iA(k)) > 0.0) ; gB = (f(iB(k)) > 0.0)
    if ((gA .eqv. gC) .and. (gB .eqv. gC)) then
      ! Uncut triangle: interior 3-pt rule on (C, A, B); QPs (2,3) are a reflection orbit.
      bA(:) = 0.0 ; bA(iA(k)) = 1.0
      bB(:) = 0.0 ; bB(iB(k)) = 1.0
      beta(:,1,k) = (SEP2_W23 * bC(:)) + ((SEP2_W16 * bA(:)) + (SEP2_W16 * bB(:)))
      beta(:,2,k) = (SEP2_W16 * bC(:)) + ((SEP2_W23 * bA(:)) + (SEP2_W16 * bB(:)))
      beta(:,3,k) = (SEP2_W16 * bC(:)) + ((SEP2_W16 * bA(:)) + (SEP2_W23 * bB(:)))
      wref(1:3,k) = SEP2_TRI3
      qp_grounded(1:3,k) = gC
      nqp(k) = 3
    elseif (gA .eqv. gB) then
      ! Center vertex separated; cyclic role binding (C, A, B).
      bA(:) = 0.0 ; bA(iA(k)) = 1.0
      bB(:) = 0.0 ; bB(iB(k)) = 1.0
      call sep2_cut_tri(fC, f(iA(k)), f(iB(k)), bC, bA, bB, gC, beta(:,:,k), wref(:,k), qp_grounded(:,k))
      nqp(k) = 7
    elseif (gB .eqv. gC) then
      ! Corner A separated; cyclic role binding (A, B, C).
      bA(:) = 0.0 ; bA(iA(k)) = 1.0
      bB(:) = 0.0 ; bB(iB(k)) = 1.0
      call sep2_cut_tri(f(iA(k)), f(iB(k)), fC, bA, bB, bC, gA, beta(:,:,k), wref(:,k), qp_grounded(:,k))
      nqp(k) = 7
    else
      ! Corner B separated; cyclic role binding (B, C, A).
      bA(:) = 0.0 ; bA(iA(k)) = 1.0
      bB(:) = 0.0 ; bB(iB(k)) = 1.0
      call sep2_cut_tri(f(iB(k)), fC, f(iA(k)), bB, bC, bA, gB, beta(:,:,k), wref(:,k), qp_grounded(:,k))
      nqp(k) = 7
    endif
  enddo

end subroutine sep2_cell_qps

!> Quadrature of one cut parent triangle: minority vertex X separated from (Y, Z) by the
!! straight zero contour of the linear deficit. QPs 1-3 sample the X-side triangle piece
!! (interior 3-pt rule, X-heavy first; 2 and 3 are a reflection orbit); QPs 4-7 sample the
!! (Y,Z)-side quad piece (2x2 tensor rule on the bilinear sub-map; pairs (4,5) and (6,7)
!! are reflection orbits). A cut through a vertex collapses the quad exactly (zero-area
!! side), so degenerate configurations need no special case.
subroutine sep2_cut_tri(fX, fY, fZ, bX, bY, bZ, gX, betaT, wrefT, gT)
  real,               intent(in)  :: fX     !< Deficit at the minority vertex [Z ~> m]
  real,               intent(in)  :: fY     !< Deficit at the first (CCW) majority vertex [Z ~> m]
  real,               intent(in)  :: fZ     !< Deficit at the second majority vertex [Z ~> m]
  real, dimension(4), intent(in)  :: bX     !< Corner-basis weights of vertex X [nondim]
  real, dimension(4), intent(in)  :: bY     !< Corner-basis weights of vertex Y [nondim]
  real, dimension(4), intent(in)  :: bZ     !< Corner-basis weights of vertex Z [nondim]
  logical,            intent(in)  :: gX     !< Grounded state of the minority vertex
  real, dimension(4,7), intent(out) :: betaT !< Corner-basis weights per (corner, QP) [nondim]
  real, dimension(7),   intent(out) :: wrefT !< Reference-space measure per QP [nondim]
  logical, dimension(7), intent(out) :: gT   !< Grounded state per QP

  real :: cY1, cX1  ! Crossing weights on edge X-Y: v1 = cX1*X + cY1*Y [nondim]
  real :: cZ2, cX2  ! Crossing weights on edge Z-X: v4 = cX2*X + cZ2*Z [nondim]
  real, dimension(4) :: b1, b4 ! Corner-basis weights of the crossings [nondim]
  real :: t1, t2, t3, t4 ! Tensor-product factors at a quad QP [nondim]
  real :: wtri      ! Per-QP measure of the triangle piece [nondim]
  integer :: k, ir, is

  ! Exact edge crossings; both complements have their own role-anchored formula so a
  ! (Y,Z) swap permutes them bitwise (denominators are exact negations of each other).
  cY1 = fX / (fX - fY) ; cX1 = fY / (fY - fX)
  cZ2 = fX / (fX - fZ) ; cX2 = fZ / (fZ - fX)
  b1(:) = (cX1 * bX(:)) + (cY1 * bY(:))
  b4(:) = (cX2 * bX(:)) + (cZ2 * bZ(:))

  ! Triangle piece (X, v1, v4); ref area = cY1*cZ2 * (1/4), the parent-triangle area.
  wtri = (cY1 * cZ2) * SEP2_TRI3
  betaT(:,1) = (SEP2_W23 * bX(:)) + ((SEP2_W16 * b1(:)) + (SEP2_W16 * b4(:)))
  betaT(:,2) = (SEP2_W16 * bX(:)) + ((SEP2_W23 * b1(:)) + (SEP2_W16 * b4(:)))
  betaT(:,3) = (SEP2_W16 * bX(:)) + ((SEP2_W16 * b1(:)) + (SEP2_W23 * b4(:)))
  wrefT(1:3) = wtri
  gT(1:3) = gX

  ! Quad piece (v1, Y, Z, v4) on the bilinear sub-map Q(r,s): r along v1->Y and v4->Z,
  ! s along v1->v4. Sub-map Jacobian (linear in r and s, derived analytically so a (Y,Z)
  ! swap maps it to J(r,1-s) bitwise; 0.5 = cross(Y-X, Z-X), twice the parent-tri area):
  !   J_sub = 0.5 * [ (1-s)*cX1*((1-r)*cZ2 + r) + s*cX2*((1-r)*cY1 + r) ]
  ! QP weight = (1/4 Gauss) * J_sub; 0.125 = 0.25 * 0.5.
  k = 3
  do ir=1,2 ; do is=1,2
    k = k + 1
    t1 = SEP2_GC(ir) * SEP2_GC(is) ; t2 = SEP2_GP(ir) * SEP2_GC(is)
    t3 = SEP2_GP(ir) * SEP2_GP(is) ; t4 = SEP2_GC(ir) * SEP2_GP(is)
    betaT(:,k) = ((t1 * b1(:)) + (t4 * b4(:))) + ((t2 * bY(:)) + (t3 * bZ(:)))
    wrefT(k) = 0.125 * ( ((SEP2_GC(is) * cX1) * ((SEP2_GC(ir) * cZ2) + SEP2_GP(ir))) + &
                         ((SEP2_GP(is) * cX2) * ((SEP2_GC(ir) * cY1) + SEP2_GP(ir))) )
    gT(k) = .not. gX
  enddo ; enddo

end subroutine sep2_cut_tri

!> Reference-space gradient of the P1 (linear) interpolant of a corner field on each of the four
!! parent triangles of the sep2_cell_qps fan. Triangle t is (C, A, B) with C the cell center, whose
!! value is the corner mean -- formed here with the same diagonal-pair grouping sep2_cell_qps uses,
!! so it is bitwise the value the partition was built on. The gradient is constant on each triangle.
!!
!! Every field that enters a quantity which branches on the SEP2 flotation state must be
!! differentiated with THIS operator, not with the bilinear (Q1) gradient. The two agree at the
!! corners, at the center, and along the cell edges, but differ inside the cell by eta times the
!! twist term (f_SW + f_NE) - (f_SE + f_NW): the bilinear x-gradient blends the north and south edge
!! differences, while the P1 gradient on the southern triangle freezes the south edge value.
!! Mixing them breaks the surface reconstruction in two ways. It stops the grounded branch
!! telescoping -- (1-r)*B[h] + P[r*h - b] is a blend (1-r)*B + r*P applied to h against a pure P
!! applied to b, which is no consistent surface at all -- and it moves the kink off the cut, since
!! S_grounded - S_floating = r*h - b identically, so a branch taken on the P1 contour while the
!! slopes come from the bilinear field makes S jump across the cut by the bilinear r*h - b evaluated
!! there. Using P1 throughout gives grad S = P[h] - P[b] grounded and (1-r)*P[h] floating, both
!! exact, with the jump identically zero. It also matches the friction, which evaluates fields at a
!! QP as sum(beta*f) -- the beta are barycentric coordinates on the parent triangle, so that sum is
!! already the P1 interpolant.
subroutine sep2_fan_gradient(f, gxi, geta)
  real, dimension(4), intent(in)  :: f    !< Corner values, ordered SW, SE, NW, NE [Z ~> m]
  real, dimension(4), intent(out) :: gxi  !< d/dxi of the P1 interpolant per triangle, ordered
                                          !! S, E, N, W as in sep2_cell_qps [Z ~> m]
  real, dimension(4), intent(out) :: geta !< d/deta of the P1 interpolant per triangle [Z ~> m]

  ! Base-corner ids (A,B) and reference coordinates of triangles S, E, N, W; these must match
  ! sep2_cell_qps exactly or the gradient will not belong to the interpolant that was cut.
  integer, dimension(4), parameter :: iA = (/ 1, 2, 4, 3 /) ! First (CCW) base corner per triangle
  integer, dimension(4), parameter :: iB = (/ 2, 4, 3, 1 /) ! Second (CCW) base corner per triangle
  real, dimension(4), parameter :: xref = (/ 0.0, 1.0, 0.0, 1.0 /) ! Corner xi [nondim]
  real, dimension(4), parameter :: yref = (/ 0.0, 0.0, 1.0, 1.0 /) ! Corner eta [nondim]
  real :: fC          ! Field value at the cell center, the corner mean [Z ~> m]
  real :: abx, aby    ! Reference-space edge vector A->B [nondim]
  real :: acx, acy    ! Reference-space edge vector A->C [nondim]
  real :: u, v        ! Field increments along A->B and A->C [Z ~> m]
  real :: Idet        ! Inverse reference-space determinant of (AB, AC) [nondim]
  integer :: t

  fC = 0.5 * ((0.5 * (f(1) + f(4))) + (0.5 * (f(2) + f(3))))
  do t=1,4
    abx = xref(iB(t)) - xref(iA(t)) ; aby = yref(iB(t)) - yref(iA(t))
    acx = 0.5 - xref(iA(t))         ; acy = 0.5 - yref(iA(t))
    u = f(iB(t)) - f(iA(t))         ; v = fC - f(iA(t))
    ! Cramer solve of g.AB = u, g.AC = v for the constant gradient g on this triangle.
    Idet = 1.0 / ((abx * acy) - (acx * aby))
    gxi(t)  = ((u * acy) - (v * aby)) * Idet
    geta(t) = ((v * abx) - (u * acx)) * Idet
  enddo

end subroutine sep2_fan_gradient

!> SEP2 subgrid basal traction for a CG action: Picard and Newton friction integrated over
!! the grounded sub-elements of the sep2_cell_qps partition. QPs inherit their piece's
!! flotation state; floating QPs contribute nothing.
subroutine CG_action_sep2_basal(CS, G, US, H, U_curr, V_curr, U_delta, V_delta, &
                                bathyT, dens_ratio, i_elem, j_elem, fB_e, use_newton, &
                                Ucontr, Vcontr, dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT, &
                                use_DG, h_nodal_cell, bed_corners, fls_cell)
  type(ice_shelf_dyn_CS), intent(in) :: CS      !< Ice shelf control structure
  type(ocean_grid_type),  intent(in) :: G       !< The grid structure
  type(unit_scale_type),  intent(in) :: US      !< Unit conversion factors
  real, dimension(2,2),   intent(in) :: H       !< Ice thickness at element corners [Z ~> m]
  real, dimension(2,2),   intent(in) :: U_curr  !< Frozen u^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_curr  !< Frozen v^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: U_delta !< Search direction du at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_delta !< Search direction dv at element corners [L T-1 ~> m s-1]
  real,                   intent(in) :: bathyT  !< Ocean bathymetry depth at tracer point [Z ~> m]
  real,                   intent(in) :: dens_ratio !< Ice density / water density [nondim]
  integer,                intent(in) :: i_elem  !< Tracer-grid i-index of the element
  integer,                intent(in) :: j_elem  !< Tracer-grid j-index of the element
  real,                   intent(in) :: fB_e    !< Element Coulomb parameter fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  logical,                intent(in) :: use_newton !< If true, include Newton basal drag correction
  real, dimension(2,2),   intent(out) :: Ucontr !< Nodal u-contributions with friction applied [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2),   intent(out) :: Vcontr !< Nodal v-contributions with friction applied [R L3 Z T-2 ~> kg m s-2]
  real,                   intent(in) :: dxCv_S  !< The cell width at the southern (v-point) edge [L ~> m]
  real,                   intent(in) :: dxCv_N  !< The cell width at the northern (v-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_W  !< The cell height at the western (u-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_E  !< The cell height at the eastern (u-point) edge [L ~> m]
  real,                   intent(in) :: IareaT  !< The inverse of the cell area at the tracer point [L-2 ~> m-2]
  logical,       optional, intent(in) :: use_DG !< If true, use DG nodal thickness and bed corners
  real, dimension(2,2), optional, intent(in) :: h_nodal_cell !< Q1 thickness at the 4 cell corners [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: bed_corners  !< Bed elevation at element corners [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: fls_cell !< Flotation deficit r*h - bed at the 4 cell
                                              !! corners (I2N_FRICTION). When present the
                                              !! partition and the effective pressure are both taken
                                              !! from this field and h_nodal_cell supplies the
                                              !! matching corner thickness [Z ~> m]

  real, dimension(4)     :: hc, bedc  ! Corner thickness and bed, flattened SW,SE,NW,NE [Z ~> m]
  real, dimension(4)     :: uc, vc    ! Corner frozen velocities [L T-1 ~> m s-1]
  real, dimension(4)     :: duc, dvc  ! Corner search directions [L T-1 ~> m s-1]
  real, dimension(4)     :: fls       ! Corner flotation deficit r*h - bed [Z ~> m]
  integer, dimension(4)  :: nqp       ! QPs per parent triangle
  real, dimension(4,7,4) :: beta      ! Corner-basis weights per (corner, QP, triangle) [nondim]
  real, dimension(7,4)   :: wref      ! Reference measure per (QP, triangle) [nondim]
  logical, dimension(7,4) :: qpg      ! Grounded state per (QP, triangle)
  real, dimension(4,7)   :: valu, valv ! Per-QP nodal contributions [R L3 Z T-2 ~> kg m s-2]
  real, dimension(4,4)   :: pu, pv    ! Per-(corner, triangle) partial sums [R L3 Z T-2 ~> kg m s-2]
  real :: b1, b2, b3, b4    ! Corner-basis weights at the QP [nondim]
  real :: mS, mN, mW, mE    ! Marginal sums: interpolation weights of the 4 cell edges [nondim]
  real :: a, d              ! Interpolated cell-edge spacings at the QP [L ~> m]
  real :: jac               ! Quadrature weight wref * (a*d) * IareaT [nondim]
  real :: hloc              ! Ice thickness at the QP [Z ~> m]
  real :: bed_sub           ! Bed elevation at the QP [Z ~> m]
  real :: u_curr_loc, v_curr_loc   ! Frozen velocity at the QP [L T-1 ~> m s-1]
  real :: u_delta_loc, v_delta_loc ! Search direction at the QP [L T-1 ~> m s-1]
  real :: unorm2_loc        ! Regularized |u^k|^2 at the QP [L2 T-2 ~> m2 s-2]
  real :: basal_coef_loc    ! Picard friction coefficient at the QP [R L2 Z T-1 ~> kg s-1]
  real :: drag_newt_loc     ! Newton drag coefficient at the QP [R Z T ~> kg m-2 s]
  real :: inner_dot_loc     ! u^k . du inner product at the QP [L2 T-2 ~> m2 s-2]
  real :: contrib           ! Per-corner quadrature contribution [nondim]
  real :: coef_prefactor    ! Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real :: min_trac_area     ! Minimum area-integrated traction floor [R L2 Z T-1 ~> kg s-1]
  real :: eps_vel2          ! Velocity regularization squared [L2 T-2 ~> m2 s-2]
  real :: fB_local          ! Coulomb fB at the QP (DG mode) [(T L-1)^CF_PostPeak]
  real :: rho_oi_ratio      ! density_ocean / density_ice [nondim]
  real :: rho_ice_g_LtoZ    ! US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  real :: rho_ocean_g_LtoZ  ! US%L_to_Z * density_ocean_avg * g_Earth [R L Z-1 T-2]
  real :: fls_loc           ! Flotation deficit at the QP [Z ~> m]
  logical :: do_DG          ! Local flag for DG mode
  logical :: do_fvsub       ! Local flag for the FV sub-element mode (fls_cell supplied)
  integer :: t, k, c

  coef_prefactor = CS%coef_prefactor(i_elem,j_elem)
  min_trac_area  = CS%min_basal_traction * G%areaT(i_elem,j_elem)
  eps_vel2 = CS%eps_glen_min**2 * ((G%dxT(i_elem,j_elem)**2) + (G%dyT(i_elem,j_elem)**2))

  do_DG = .false.
  if (present(use_DG)) do_DG = use_DG
  do_fvsub = present(fls_cell)
  if (do_fvsub) then
    ! FV sub-element mode: the partition, the sub-element thickness and the effective pressure all
    ! come from the two corner fields built by build_corner_flotation_fields, so the grounding line
    ! the friction sees is the same one the driving stress sees. The bed is implicit in fls.
    rho_ocean_g_LtoZ = US%L_to_Z * CS%density_ocean_avg * CS%g_Earth
    hc(1) = h_nodal_cell(1,1) ; hc(2) = h_nodal_cell(2,1)
    hc(3) = h_nodal_cell(1,2) ; hc(4) = h_nodal_cell(2,2)
    bedc(:) = bathyT
  elseif (do_DG) then
    rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
    rho_ice_g_LtoZ = US%L_to_Z * CS%density_ice * CS%g_Earth
    hc(1) = h_nodal_cell(1,1) ; hc(2) = h_nodal_cell(2,1)
    hc(3) = h_nodal_cell(1,2) ; hc(4) = h_nodal_cell(2,2)
    bedc(1) = bed_corners(1,1) ; bedc(2) = bed_corners(2,1)
    bedc(3) = bed_corners(1,2) ; bedc(4) = bed_corners(2,2)
  else
    hc(1) = H(1,1) ; hc(2) = H(2,1) ; hc(3) = H(1,2) ; hc(4) = H(2,2)
    bedc(:) = bathyT
  endif
  uc(1)  = U_curr(1,1)  ; uc(2)  = U_curr(2,1)  ; uc(3)  = U_curr(1,2)  ; uc(4)  = U_curr(2,2)
  vc(1)  = V_curr(1,1)  ; vc(2)  = V_curr(2,1)  ; vc(3)  = V_curr(1,2)  ; vc(4)  = V_curr(2,2)
  duc(1) = U_delta(1,1) ; duc(2) = U_delta(2,1) ; duc(3) = U_delta(1,2) ; duc(4) = U_delta(2,2)
  dvc(1) = V_delta(1,1) ; dvc(2) = V_delta(2,1) ; dvc(3) = V_delta(1,2) ; dvc(4) = V_delta(2,2)

  if (do_fvsub) then
    ! The corner deficit is already consistent with hc (both carried by one weight set over one cell
    ! set), and MIN_H_SHELF was applied at the cell centers, so no clamp is reapplied here.
    fls(1) = fls_cell(1,1) ; fls(2) = fls_cell(2,1)
    fls(3) = fls_cell(1,2) ; fls(4) = fls_cell(2,2)
  else
    ! Unclamped bilinear deficit defines the partition; magnitudes keep the min_h clamp below.
    fls(:) = (dens_ratio * hc(:)) - bedc(:)
  endif
  call sep2_cell_qps(fls, nqp, beta, wref, qpg)

  do t=1,4
    valu(:,1:nqp(t)) = 0.0 ; valv(:,1:nqp(t)) = 0.0
    do k=1,nqp(t)
      if (.not. qpg(k,t)) cycle ! floating piece: no basal traction
      b1 = beta(1,k,t) ; b2 = beta(2,k,t) ; b3 = beta(3,k,t) ; b4 = beta(4,k,t)
      mS = b1 + b2 ; mN = b3 + b4 ; mW = b1 + b3 ; mE = b2 + b4
      a = (dxCv_S * mS) + (dxCv_N * mN)
      d = (dyCu_W * mW) + (dyCu_E * mE)
      jac = (wref(k,t) * (a * d)) * IareaT

      hloc = ((b1 * hc(1)) + (b4 * hc(4))) + ((b2 * hc(2)) + (b3 * hc(3)))
      if (do_fvsub) then
        ! The corner-basis weights beta are the barycentric coordinates of the QP in its parent
        ! triangle (the center weights bC are 1/4 each, so they reproduce the cell-center value), so
        ! this sum is the P1-on-the-fan interpolant -- the very field whose zero contour sep2_cut_tri
        ! cut. Hence fls_loc >= 0 at every grounded QP by construction, and the effective pressure
        ! below can never see a negative argument.
        fls_loc = ((b1 * fls(1)) + (b4 * fls(4))) + ((b2 * fls(2)) + (b3 * fls(3)))
        bed_sub = bathyT
      elseif (do_DG) then
        hloc = max(hloc, CS%min_h_shelf)
        bed_sub = ((b1 * bedc(1)) + (b4 * bedc(4))) + ((b2 * bedc(2)) + (b3 * bedc(3)))
      else
        bed_sub = bathyT
      endif
      u_curr_loc  = ((b1 * uc(1))  + (b4 * uc(4)))  + ((b2 * uc(2))  + (b3 * uc(3)))
      v_curr_loc  = ((b1 * vc(1))  + (b4 * vc(4)))  + ((b2 * vc(2))  + (b3 * vc(3)))
      u_delta_loc = ((b1 * duc(1)) + (b4 * duc(4))) + ((b2 * duc(2)) + (b3 * duc(3)))
      v_delta_loc = ((b1 * dvc(1)) + (b4 * dvc(4))) + ((b2 * dvc(2)) + (b3 * dvc(3)))

      unorm2_loc = ((u_curr_loc**2) + (v_curr_loc**2)) + eps_vel2

      if (do_fvsub .and. CS%CoulombFriction) then
        fB_local = compute_fB_from_N( &
            subgrid_effective_pressure(fls_loc, hloc, dens_ratio, rho_ocean_g_LtoZ), &
            CS%C_basal_friction(i_elem,j_elem), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
            CS%CF_PostPeak, CS%n_basal_fric)
      elseif (do_DG .and. CS%CoulombFriction) then
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

      do c=1,4
        contrib = jac * beta(c,k,t)
        valu(c,k) = contrib * (basal_coef_loc * u_delta_loc)
        valv(c,k) = contrib * (basal_coef_loc * v_delta_loc)
        if (use_newton) then
          valu(c,k) = valu(c,k) + (contrib * (drag_newt_loc * u_curr_loc * inner_dot_loc))
          valv(c,k) = valv(c,k) + (contrib * (drag_newt_loc * v_curr_loc * inner_dot_loc))
        endif
      enddo
    enddo

    ! Orbit-grouped QP sums: (2,3), (4,5) and (6,7) are reflection pairs.
    if (nqp(t) == 3) then
      do c=1,4
        pu(c,t) = valu(c,1) + (valu(c,2) + valu(c,3))
        pv(c,t) = valv(c,1) + (valv(c,2) + valv(c,3))
      enddo
    else
      do c=1,4
        pu(c,t) = (valu(c,1) + (valu(c,2) + valu(c,3))) + &
                  ((valu(c,4) + valu(c,5)) + (valu(c,6) + valu(c,7)))
        pv(c,t) = (valv(c,1) + (valv(c,2) + valv(c,3))) + &
                  ((valv(c,4) + valv(c,5)) + (valv(c,6) + valv(c,7)))
      enddo
    endif
  enddo

  ! Role-grouped cross-triangle reduction: each corner takes each of the roles
  ! (A, B, farA, farB) exactly once over the 4 triangles (S=1, E=2, N=3, W=4).
  Ucontr(1,1) = (pu(1,1) + pu(1,4)) + (pu(1,2) + pu(1,3)) ! SW: (S+W)+(E+N)
  Ucontr(2,1) = (pu(2,2) + pu(2,1)) + (pu(2,3) + pu(2,4)) ! SE: (E+S)+(N+W)
  Ucontr(1,2) = (pu(3,4) + pu(3,3)) + (pu(3,1) + pu(3,2)) ! NW: (W+N)+(S+E)
  Ucontr(2,2) = (pu(4,3) + pu(4,2)) + (pu(4,4) + pu(4,1)) ! NE: (N+E)+(W+S)
  Vcontr(1,1) = (pv(1,1) + pv(1,4)) + (pv(1,2) + pv(1,3))
  Vcontr(2,1) = (pv(2,2) + pv(2,1)) + (pv(2,3) + pv(2,4))
  Vcontr(1,2) = (pv(3,4) + pv(3,3)) + (pv(3,1) + pv(3,2))
  Vcontr(2,2) = (pv(4,3) + pv(4,2)) + (pv(4,4) + pv(4,1))

end subroutine CG_action_sep2_basal

!> SEP2 subgrid basal traction for the preconditioner diagonal: same partition and
!! quadrature as CG_action_sep2_basal, with squared basis weights and per-block
!! Newton terms (u^2 for the u-block, v^2 for the v-block).
subroutine CG_diagonal_sep2_basal(CS, G, US, H, U_curr, V_curr, &
                                  bathyT, dens_ratio, i_elem, j_elem, fB_e, u_diag, v_diag, &
                                  dxCv_S, dxCv_N, dyCu_W, dyCu_E, IareaT, &
                                  use_DG, h_nodal_cell, bed_corners, fls_cell)
  type(ice_shelf_dyn_CS), intent(in) :: CS      !< Ice shelf control structure
  type(ocean_grid_type),  intent(in) :: G       !< The grid structure
  type(unit_scale_type),  intent(in) :: US      !< Unit conversion factors
  real, dimension(2,2),   intent(in) :: H       !< Ice thickness at element corners [Z ~> m]
  real, dimension(2,2),   intent(in) :: U_curr  !< Frozen u^k at element corners [L T-1 ~> m s-1]
  real, dimension(2,2),   intent(in) :: V_curr  !< Frozen v^k at element corners [L T-1 ~> m s-1]
  real,                   intent(in) :: bathyT  !< Ocean bathymetry depth at tracer point [Z ~> m]
  real,                   intent(in) :: dens_ratio !< Ice density / water density [nondim]
  integer,                intent(in) :: i_elem  !< Tracer-grid i-index of the element
  integer,                intent(in) :: j_elem  !< Tracer-grid j-index of the element
  real,                   intent(in) :: fB_e    !< Element Coulomb parameter fB; 0 for Weertman [(T L-1)^CF_PostPeak]
  real, dimension(2,2),   intent(out) :: u_diag !< Nodal u-diagonal entries [R L2 Z T-1 ~> kg s-1]
  real, dimension(2,2),   intent(out) :: v_diag !< Nodal v-diagonal entries [R L2 Z T-1 ~> kg s-1]
  real,                   intent(in) :: dxCv_S  !< The cell width at the southern (v-point) edge [L ~> m]
  real,                   intent(in) :: dxCv_N  !< The cell width at the northern (v-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_W  !< The cell height at the western (u-point) edge [L ~> m]
  real,                   intent(in) :: dyCu_E  !< The cell height at the eastern (u-point) edge [L ~> m]
  real,                   intent(in) :: IareaT  !< The inverse of the cell area at the tracer point [L-2 ~> m-2]
  logical,       optional, intent(in) :: use_DG !< If true, use DG nodal thickness and bed corners
  real, dimension(2,2), optional, intent(in) :: h_nodal_cell !< Q1 thickness at the 4 cell corners [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: bed_corners  !< Bed elevation at element corners [Z ~> m]
  real, dimension(2,2), optional, intent(in) :: fls_cell !< Flotation deficit r*h - bed at the 4 cell
                                              !! corners (I2N_FRICTION); see
                                              !! CG_action_sep2_basal [Z ~> m]

  real, dimension(4)     :: hc, bedc  ! Corner thickness and bed, flattened SW,SE,NW,NE [Z ~> m]
  real, dimension(4)     :: uc, vc    ! Corner frozen velocities [L T-1 ~> m s-1]
  real, dimension(4)     :: fls       ! Corner flotation deficit r*h - bed [Z ~> m]
  integer, dimension(4)  :: nqp       ! QPs per parent triangle
  real, dimension(4,7,4) :: beta      ! Corner-basis weights per (corner, QP, triangle) [nondim]
  real, dimension(7,4)   :: wref      ! Reference measure per (QP, triangle) [nondim]
  logical, dimension(7,4) :: qpg      ! Grounded state per (QP, triangle)
  real, dimension(4,7)   :: valu, valv ! Per-QP nodal diagonal contributions [R L2 Z T-1 ~> kg s-1]
  real, dimension(4,4)   :: pu, pv    ! Per-(corner, triangle) partial sums [R L2 Z T-1 ~> kg s-1]
  real :: b1, b2, b3, b4    ! Corner-basis weights at the QP [nondim]
  real :: mS, mN, mW, mE    ! Marginal sums: interpolation weights of the 4 cell edges [nondim]
  real :: a, d              ! Interpolated cell-edge spacings at the QP [L ~> m]
  real :: jac               ! Quadrature weight wref * (a*d) * IareaT [nondim]
  real :: hloc              ! Ice thickness at the QP [Z ~> m]
  real :: bed_sub           ! Bed elevation at the QP [Z ~> m]
  real :: u_curr_loc, v_curr_loc ! Frozen velocity at the QP [L T-1 ~> m s-1]
  real :: unorm2_loc        ! Regularized |u^k|^2 at the QP [L2 T-2 ~> m2 s-2]
  real :: basal_coef_loc    ! Picard friction coefficient at the QP [R L2 Z T-1 ~> kg s-1]
  real :: drag_newt_loc     ! Newton drag coefficient at the QP [R Z T ~> kg m-2 s]
  real :: contrib           ! Per-corner quadrature contribution [nondim]
  real :: coef_prefactor    ! Pre-computed area * C_basal_friction * L_T_to_m_s [R L2 Z T-1 ~> kg s-1]
  real :: min_trac_area     ! Minimum area-integrated traction floor [R L2 Z T-1 ~> kg s-1]
  real :: eps_vel2          ! Velocity regularization squared [L2 T-2 ~> m2 s-2]
  real :: fB_local          ! Coulomb fB at the QP (DG mode) [(T L-1)^CF_PostPeak]
  real :: rho_oi_ratio      ! density_ocean / density_ice [nondim]
  real :: rho_ice_g_LtoZ    ! US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  real :: rho_ocean_g_LtoZ  ! US%L_to_Z * density_ocean_avg * g_Earth [R L Z-1 T-2]
  real :: fls_loc           ! Flotation deficit at the QP [Z ~> m]
  logical :: do_DG          ! Local flag for DG mode
  logical :: do_fvsub       ! Local flag for the FV sub-element mode (fls_cell supplied)
  integer :: t, k, c

  coef_prefactor = CS%coef_prefactor(i_elem,j_elem)
  min_trac_area  = CS%min_basal_traction * G%areaT(i_elem,j_elem)
  eps_vel2 = CS%eps_glen_min**2 * ((G%dxT(i_elem,j_elem)**2) + (G%dyT(i_elem,j_elem)**2))

  do_DG = .false.
  if (present(use_DG)) do_DG = use_DG
  do_fvsub = present(fls_cell)
  if (do_fvsub) then
    ! FV sub-element mode: the partition, the sub-element thickness and the effective pressure all
    ! come from the two corner fields built by build_corner_flotation_fields, so the grounding line
    ! the friction sees is the same one the driving stress sees. The bed is implicit in fls.
    rho_ocean_g_LtoZ = US%L_to_Z * CS%density_ocean_avg * CS%g_Earth
    hc(1) = h_nodal_cell(1,1) ; hc(2) = h_nodal_cell(2,1)
    hc(3) = h_nodal_cell(1,2) ; hc(4) = h_nodal_cell(2,2)
    bedc(:) = bathyT
  elseif (do_DG) then
    rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
    rho_ice_g_LtoZ = US%L_to_Z * CS%density_ice * CS%g_Earth
    hc(1) = h_nodal_cell(1,1) ; hc(2) = h_nodal_cell(2,1)
    hc(3) = h_nodal_cell(1,2) ; hc(4) = h_nodal_cell(2,2)
    bedc(1) = bed_corners(1,1) ; bedc(2) = bed_corners(2,1)
    bedc(3) = bed_corners(1,2) ; bedc(4) = bed_corners(2,2)
  else
    hc(1) = H(1,1) ; hc(2) = H(2,1) ; hc(3) = H(1,2) ; hc(4) = H(2,2)
    bedc(:) = bathyT
  endif
  uc(1) = U_curr(1,1) ; uc(2) = U_curr(2,1) ; uc(3) = U_curr(1,2) ; uc(4) = U_curr(2,2)
  vc(1) = V_curr(1,1) ; vc(2) = V_curr(2,1) ; vc(3) = V_curr(1,2) ; vc(4) = V_curr(2,2)

  fls(:) = (dens_ratio * hc(:)) - bedc(:)
  call sep2_cell_qps(fls, nqp, beta, wref, qpg)

  do t=1,4
    valu(:,1:nqp(t)) = 0.0 ; valv(:,1:nqp(t)) = 0.0
    do k=1,nqp(t)
      if (.not. qpg(k,t)) cycle ! floating piece: no basal traction
      b1 = beta(1,k,t) ; b2 = beta(2,k,t) ; b3 = beta(3,k,t) ; b4 = beta(4,k,t)
      mS = b1 + b2 ; mN = b3 + b4 ; mW = b1 + b3 ; mE = b2 + b4
      a = (dxCv_S * mS) + (dxCv_N * mN)
      d = (dyCu_W * mW) + (dyCu_E * mE)
      jac = (wref(k,t) * (a * d)) * IareaT

      hloc = ((b1 * hc(1)) + (b4 * hc(4))) + ((b2 * hc(2)) + (b3 * hc(3)))
      if (do_fvsub) then
        ! The corner-basis weights beta are the barycentric coordinates of the QP in its parent
        ! triangle (the center weights bC are 1/4 each, so they reproduce the cell-center value), so
        ! this sum is the P1-on-the-fan interpolant -- the very field whose zero contour sep2_cut_tri
        ! cut. Hence fls_loc >= 0 at every grounded QP by construction, and the effective pressure
        ! below can never see a negative argument.
        fls_loc = ((b1 * fls(1)) + (b4 * fls(4))) + ((b2 * fls(2)) + (b3 * fls(3)))
        bed_sub = bathyT
      elseif (do_DG) then
        hloc = max(hloc, CS%min_h_shelf)
        bed_sub = ((b1 * bedc(1)) + (b4 * bedc(4))) + ((b2 * bedc(2)) + (b3 * bedc(3)))
      else
        bed_sub = bathyT
      endif
      u_curr_loc = ((b1 * uc(1)) + (b4 * uc(4))) + ((b2 * uc(2)) + (b3 * uc(3)))
      v_curr_loc = ((b1 * vc(1)) + (b4 * vc(4))) + ((b2 * vc(2)) + (b3 * vc(3)))

      unorm2_loc = ((u_curr_loc**2) + (v_curr_loc**2)) + eps_vel2

      if (do_fvsub .and. CS%CoulombFriction) then
        fB_local = compute_fB_from_N( &
            subgrid_effective_pressure(fls_loc, hloc, dens_ratio, rho_ocean_g_LtoZ), &
            CS%C_basal_friction(i_elem,j_elem), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
            CS%CF_PostPeak, CS%n_basal_fric)
      elseif (do_DG .and. CS%CoulombFriction) then
        fB_local = compute_fB_local(hloc, bed_sub, rho_oi_ratio, rho_ice_g_LtoZ, &
            CS%C_basal_friction(i_elem,j_elem), CS%alpha_coulomb, CS%CF_Max, CS%CF_MinN, &
            CS%CF_PostPeak, CS%n_basal_fric)
      else
        fB_local = fB_e
      endif

      call compute_basal_coef(unorm2_loc, coef_prefactor, min_trac_area, fB_local, &
          CS%n_basal_fric, CS%CoulombFriction, CS%CF_PostPeak, US%L_T_to_m_s, .true., &
          basal_coef_loc, drag_newt_loc)

      do c=1,4
        contrib = jac * (beta(c,k,t)**2)
        if (CS%doing_newton) then
          valu(c,k) = contrib * (basal_coef_loc + (drag_newt_loc * u_curr_loc**2))
          valv(c,k) = contrib * (basal_coef_loc + (drag_newt_loc * v_curr_loc**2))
        else
          valu(c,k) = contrib * basal_coef_loc
          valv(c,k) = contrib * basal_coef_loc
        endif
      enddo
    enddo

    ! Orbit-grouped QP sums: (2,3), (4,5) and (6,7) are reflection pairs.
    if (nqp(t) == 3) then
      do c=1,4
        pu(c,t) = valu(c,1) + (valu(c,2) + valu(c,3))
        pv(c,t) = valv(c,1) + (valv(c,2) + valv(c,3))
      enddo
    else
      do c=1,4
        pu(c,t) = (valu(c,1) + (valu(c,2) + valu(c,3))) + &
                  ((valu(c,4) + valu(c,5)) + (valu(c,6) + valu(c,7)))
        pv(c,t) = (valv(c,1) + (valv(c,2) + valv(c,3))) + &
                  ((valv(c,4) + valv(c,5)) + (valv(c,6) + valv(c,7)))
      enddo
    endif
  enddo

  ! Role-grouped cross-triangle reduction (see CG_action_sep2_basal).
  u_diag(1,1) = (pu(1,1) + pu(1,4)) + (pu(1,2) + pu(1,3)) ! SW: (S+W)+(E+N)
  u_diag(2,1) = (pu(2,2) + pu(2,1)) + (pu(2,3) + pu(2,4)) ! SE: (E+S)+(N+W)
  u_diag(1,2) = (pu(3,4) + pu(3,3)) + (pu(3,1) + pu(3,2)) ! NW: (W+N)+(S+E)
  u_diag(2,2) = (pu(4,3) + pu(4,2)) + (pu(4,4) + pu(4,1)) ! NE: (N+E)+(W+S)
  v_diag(1,1) = (pv(1,1) + pv(1,4)) + (pv(1,2) + pv(1,3))
  v_diag(2,1) = (pv(2,2) + pv(2,1)) + (pv(2,3) + pv(2,4))
  v_diag(1,2) = (pv(3,4) + pv(3,3)) + (pv(3,1) + pv(3,2))
  v_diag(2,2) = (pv(4,3) + pv(4,2)) + (pv(4,4) + pv(4,1))

end subroutine CG_diagonal_sep2_basal

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
      call interpolate_H_to_B_DG(G, CS%h_nodal, ISS%hmask, &
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
      ! This runs over every ice cell, floating ones included, where h < Hf makes the unfloored
      ! effective pressure negative. compute_fB_from_N clamps it at CF_MinN and returns
      ! FB_NO_COULOMB_DRAG if that leaves it at zero, so CF_MinN = 0 is safe here.
      Hf = max((CS%density_ocean_avg/CS%density_ice) * CS%bed_elev(i,j), 0.0)
      fN = US%L_to_Z*(CS%density_ice * CS%g_Earth) * (max(ISS%h_shelf(i,j), CS%min_h_shelf) - Hf)
      CS%fB_elem(i,j) = compute_fB_from_N(fN, CS%C_basal_friction(i,j), CS%alpha_coulomb, &
          CS%CF_Max, CS%CF_MinN, CS%CF_PostPeak, CS%n_basal_fric)
    else
      CS%fB_elem(i,j) = 0.0
    endif
  enddo ; enddo

end subroutine calc_shelf_basal_prefactors

!> Pre-compute the nodal basal-friction prefactors for CISM_FRICTION='local'. C_basal_friction is a
!! static bed property defined under grounded, floating, and ice-free cells alike, so the nodal C is
!! an area-weighted average over all four cells around the node -- the nodal C does not change as the
!! (grounded) ice front moves across the node. The Coulomb effective pressure is instead an ice-state
!! quantity: with CISM_NODAL_EFFECPRESS it is formed and capped in each cell and then averaged over
!! all four in-domain cells (ice-free cells contributing zero), as CISM does; otherwise the thickness
!! and bed are averaged over the ice-covered cells alone and N is formed from those means. These are
!! velocity-independent, so they are computed once per solve alongside calc_shelf_basal_prefactors.
!! Node (I,J) is the NE corner of cell (i,j); its four surrounding cells are
!! (I,J),(I+1,J),(I,J+1),(I+1,J+1).
subroutine calc_shelf_basal_prefactors_node(CS, ISS, G, US)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< Ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< Ice shelf state (hmask, h_shelf)
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure
  type(unit_scale_type),  intent(in)    :: US  !< Unit conversion factors

  real :: rho_oi_ratio   ! density_ocean_avg / density_ice [nondim]
  real :: rho_ice_g_LtoZ ! US%L_to_Z * density_ice * g_Earth [R L Z-1 T-2]
  real :: asum_all       ! Sum of all four in-domain cell areas at the node [L2 ~> m2]
  real :: asum_ice       ! Sum of the ice-covered cell areas at the node [L2 ~> m2]
  real :: w              ! Area weight of one cell [L2 ~> m2]
  real, dimension(4) :: aq ! Per-cell area weight at the node, by corner slot [L2 ~> m2]
  real, dimension(4) :: cq ! Per-cell area-weighted C_basal_friction [R L Z T-2 (s m-1)^n L2]
  real, dimension(4) :: iq ! Per-cell area weight where ice is present [L2 ~> m2]
  real, dimension(4) :: nq ! Per-cell area-weighted effective pressure [R Z L T-2 L2]
  integer :: kq            ! Corner slot index: 1=SW, 2=SE, 3=NW, 4=NE
  real :: Cw             ! Area-weighted sum of C_basal_friction over all four cells [R L Z T-2 (s m-1)^n L2]
  real :: Nw             ! Area-weighted sum of the cell effective pressures [R Z L T-2 L2]
  real :: hw             ! Area-weighted sum of ice thickness over ice cells [Z L2 ~> m3]
  real :: bw             ! Area-weighted sum of bed elevation over ice cells [Z L2 ~> m3]
  real :: C_n            ! Nodal area-weighted C_basal_friction [R L Z T-2 (s m-1)^n]
  real :: N_n            ! Nodal area-weighted effective pressure [R Z L T-2 ~> Pa]
  real :: h_n            ! Nodal area-weighted ice thickness [Z ~> m]
  real :: bed_n          ! Nodal area-weighted bed elevation [Z ~> m]
  logical :: ice_here    ! True if this cell holds ice
  integer :: i, j, ii, jj, ic, jc
  integer :: i_off, j_off, gisc, gjsc, giec, gjec

  rho_oi_ratio   = CS%density_ocean_avg / CS%density_ice
  rho_ice_g_LtoZ = US%L_to_Z * (CS%density_ice * CS%g_Earth)
  i_off = G%idg_offset ; j_off = G%jdg_offset
  gisc = 1 ; gjsc = 1 ; giec = G%domain%niglobal ; gjec = G%domain%njglobal

  do J=G%jsd,G%jed-1 ; do I=G%isd,G%ied-1
    asum_all = 0.0 ; asum_ice = 0.0 ; Cw = 0.0 ; Nw = 0.0 ; hw = 0.0 ; bw = 0.0
    C_n = 0.0 ; N_n = 0.0
    aq(:) = 0.0 ; cq(:) = 0.0 ; iq(:) = 0.0 ; nq(:) = 0.0
    do jj=0,1 ; jc = J+jj ; do ii=0,1 ; ic = I+ii
      kq = 1 + ii + 2*jj   ! 1=SW, 2=SE, 3=NW, 4=NE
      ! Skip across-wall halo cells: beyond a non-reentrant wall C_basal_friction is not read from
      ! file and pass_var does not fill it, so it keeps the (large) allocate default and would poison
      ! the nodal C average at wall nodes, inflating the wall-node drag and breaking meridional
      ! symmetry (MISMIP3D y-velocity/y-slope). In-domain ice-free land cells are kept (valid C).
      if (.not. CS%reentrant_x) then
        if ((ic+i_off < gisc) .or. (ic+i_off > giec)) cycle
      endif
      if (.not. CS%reentrant_y) then
        if ((jc+j_off < gjsc) .or. (jc+j_off > gjec)) cycle
      endif
      w = G%areaT(ic,jc)
      ice_here = (ISS%hmask(ic,jc) == 1 .or. ISS%hmask(ic,jc) == 3)
      ! C is a static bed property under floating and ice-free ice alike: average over all in-domain
      ! cells so the nodal C is fixed by geometry, not by where the ice front happens to be.
      aq(kq) = w
      cq(kq) = w*CS%C_basal_friction(ic,jc)
      if (ice_here) iq(kq) = w
      ! CISM order of operations (glissade_basal_traction, calc_effective_pressure): form N in the
      ! cell and cap it to [0, overburden] there, then stagger N itself to the node over ALL four
      ! in-domain cells. An ice-free cell has no overburden and so contributes N = 0, exactly as
      ! CISM's glissade_stagger with stagger_margin = 0 does. Including the ice-free cells is what
      ! makes the nodal N continuous as a cell gains or loses thin ice, and capping per cell stops a
      ! deeply floating neighbor from pulling the nodal average below zero without bound.
      if (ice_here) &
        nq(kq) = w*coulomb_effective_pressure(max(ISS%h_shelf(ic,jc), CS%min_h_shelf), &
                      CS%bed_elev(ic,jc), rho_oi_ratio, rho_ice_g_LtoZ, 0.0)
    enddo ; enddo
    ! Reduce over the four cells by opposite pairs, SW+NE then SE+NW. A quarter turn
    ! permutes the four cyclically, so this pairing is bit-for-bit rotation-invariant
    ! while a running accumulation over the loop order is not. Only Nw actually varies
    ! from cell to cell here under a uniform grid and a constant C, but all four are
    ! grouped the same way so the routine stays invariant on a stretched grid too.
    asum_all = (aq(1) + aq(4)) + (aq(2) + aq(3))
    Cw       = (cq(1) + cq(4)) + (cq(2) + cq(3))
    asum_ice = (iq(1) + iq(4)) + (iq(2) + iq(3))
    Nw       = (nq(1) + nq(4)) + (nq(2) + nq(3))
    ! Nodal control volume for the local drag, as lumped corner areas (0.25*areaT per cell). With
    ! LOCAL_NODE_FULL_AREA this is the full dual-cell area of the in-domain cells, which is what CISM
    ! adds to the diagonal (dx*dy) at every active vertex; otherwise it is restricted to the
    ! ice-covered cells. Neither is areaBu: at a domain-edge node the across-wall halo cells are out
    ! of the domain, and counting them would give the boundary node a full-strength drag while its
    ! driving stress and viscous terms are integrated over only the interior (half) control volume.
    ! That mismatch slows the wall-node along-flow velocity relative to the interior, breaking the
    ! meridional symmetry of channel configs (MISMIP3D) -- spurious y-velocity and y-surface-slope
    ! near the y walls. The two choices agree wherever every in-domain cell at the node holds ice, so
    ! they differ only at ice margins, and never at the walls of an ice-filled channel.
    if (asum_ice > 0.0) then
      CS%area_node(I,J) = 0.25 * asum_all
    else
      ! No ice anywhere around the node: it is not an active vertex and carries no drag (CISM vmask).
      CS%area_node(I,J) = 0.0
    endif
    if (CS%area_node(I,J) > 0.0) then
      C_n = Cw/asum_all
      CS%coef_prefactor_node(I,J) = (CS%area_node(I,J) * C_n) * US%L_T_to_m_s
    else
      CS%coef_prefactor_node(I,J) = 0.0
    endif
    CS%fB_node(I,J) = 0.0
    if (CS%CoulombFriction .and. asum_ice > 0.0) then
      N_n = Nw/asum_all
      if (N_n > 0.0) then
        CS%fB_node(I,J) = compute_fB_from_N(N_n, C_n, CS%alpha_coulomb, CS%CF_Max, 0.0, &
            CS%CF_PostPeak, CS%n_basal_fric)
      else
        ! A vanishing effective pressure gives no Coulomb drag at all. Zero the prefactor rather
        ! than pass an infinite fB through the sliding law; MIN_BASAL_TRACTION is still applied
        ! downstream, so a partly grounded node keeps its floor under BETA_LIMIT_ABSOLUTE.
        CS%coef_prefactor_node(I,J) = 0.0
      endif
    endif
  enddo ; enddo

end subroutine calc_shelf_basal_prefactors_node

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
  real :: h_ip, bed_ip   ! Ice thickness and bed elevation at a sub-IP [Z ~> m]
  real :: bed_corners(2,2) ! Bed elevation at the 4 B-grid corners of a cell [Z ~> m]
  real :: H_corners(2,2)   ! Ice thickness at the 4 B-grid corners (non-DG path) [Z ~> m]
  real :: fls_corners(2,2) ! Flotation deficit r*h - bed at the 4 B-grid corners [Z ~> m]
  logical :: fv_sub        ! True on the FV sub-element path (CS%i2n_friction)
  real :: d_min, d_max   ! Min/max over the 4 corners of the unclamped flotation
                         ! deficit r*h - bed [Z ~> m]
  real :: bed_min        ! Min bed elevation over the 4 corners [Z ~> m]
  integer :: i, j, isub, jsub, iq, jq, n_total, n_grounded
  integer :: isc, iec, jsc, jec
  real :: g_ip                ! Flotation deficit r*h_ip - bed_ip at a sub-IP [Z ~> m]
  real, dimension(4)     :: fls_gf   ! Corner flotation deficit, flattened SW,SE,NW,NE [Z ~> m]
  integer, dimension(4)  :: nqp_gf   ! SEP2 QPs per parent triangle
  real, dimension(4,7,4) :: beta_gf  ! SEP2 corner-basis weights per (corner, QP, triangle) [nondim]
  real, dimension(7,4)   :: wref_gf  ! SEP2 reference measure per (QP, triangle) [nondim]
  logical, dimension(7,4) :: qpg_gf  ! SEP2 grounded state per (QP, triangle)
  real, dimension(7) :: vg_gf, vt_gf ! Per-QP grounded and total Jacobian weights [L2 ~> m2]
  real, dimension(4) :: pg_gf, pt_gf ! Per-triangle grounded and total weight sums [L2 ~> m2]
  real :: b1_gf, b2_gf, b3_gf, b4_gf ! Corner-basis weights at a SEP2 QP [nondim]
  real :: mS_gf, mN_gf, mW_gf, mE_gf ! Marginal edge-interpolation weights [nondim]
  real :: a_gf, d_gf                 ! Interpolated cell-edge spacings at the QP [L ~> m]
  real :: w_ground, w_total          ! Grounded and total Jacobian weights of the cell [L2 ~> m2]
  real, dimension(4,7) :: vxg_gf, vxt_gf ! Per-(corner,QP) grounded and total weights [L2 ~> m2]
  real, dimension(4,4) :: pxg_gf, pxt_gf ! Per-(corner,triangle) grounded and total sums [L2 ~> m2]
  real, dimension(4)   :: xg_gf, xt_gf   ! Per-corner grounded and total weights [L2 ~> m2]
  real, dimension(2,2) :: xi_num, xi_den ! Per-corner floating and total basis weights [nondim]
  real :: wq_gf                          ! Basis weight of a corner at a sub-IP [nondim]
  integer :: tq, kq, cq, ia, ib

  if (.not. CS%GL_regularize) return

  if (CS%use_DG_thickness) then
    endif

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  n_total = CS%n_sub_regularize * CS%n_sub_regularize * 4
  fv_sub = CS%i2n_friction

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  do j=jsc,jec ; do i=isc,iec
    if (ISS%hmask(i,j) /= 1 .and. ISS%hmask(i,j) /= 3) cycle

    ! Gather corner values for bed (and H if non-DG) at this cell, and the flotation deficit that
    ! every grounding test below keys off.
    if (CS%use_DG_thickness) then
      bed_corners(1,1) = CS%bed_node(i-1,j-1) ; bed_corners(2,1) = CS%bed_node(i,j-1)
      bed_corners(1,2) = CS%bed_node(i-1,j  ) ; bed_corners(2,2) = CS%bed_node(i,j  )
      fls_corners(1,1) = (rhoi_rhow*CS%h_nodal(i,j,1,1)) - bed_corners(1,1)
      fls_corners(2,1) = (rhoi_rhow*CS%h_nodal(i,j,2,1)) - bed_corners(2,1)
      fls_corners(1,2) = (rhoi_rhow*CS%h_nodal(i,j,1,2)) - bed_corners(1,2)
      fls_corners(2,2) = (rhoi_rhow*CS%h_nodal(i,j,2,2)) - bed_corners(2,2)
    elseif (fv_sub) then
      ! FV sub-element path: both corner fields come from build_corner_flotation_fields, which
      ! carried h and r*h-bed to the corners with one weight set, so fls = r*H - bed holds here and
      ! the bed never has to be reconstructed.
      H_corners(1,1) = CS%H_corner(i-1,j-1) ; H_corners(2,1) = CS%H_corner(i,j-1)
      H_corners(1,2) = CS%H_corner(i-1,j  ) ; H_corners(2,2) = CS%H_corner(i,j  )
      fls_corners(1,1) = CS%fls_corner(i-1,j-1) ; fls_corners(2,1) = CS%fls_corner(i,j-1)
      fls_corners(1,2) = CS%fls_corner(i-1,j  ) ; fls_corners(2,2) = CS%fls_corner(i,j  )
      bed_corners(:,:) = CS%bed_elev(i,j)  ! only reached by the min_h early-out, skipped when fv_sub
    else
      ! Non-DG: bed is cell-constant in the existing GL detection logic; mirror that
      ! by setting all 4 corner values to bed_elev(i,j).
      bed_corners(:,:) = CS%bed_elev(i,j)
      H_corners(1,1) = H_node(i-1,j-1) ; H_corners(2,1) = H_node(i,j-1)
      H_corners(1,2) = H_node(i-1,j  ) ; H_corners(2,2) = H_node(i,j  )
      fls_corners(1,1) = (rhoi_rhow*H_corners(1,1)) - bed_corners(1,1)
      fls_corners(2,1) = (rhoi_rhow*H_corners(2,1)) - bed_corners(2,1)
      fls_corners(1,2) = (rhoi_rhow*H_corners(1,2)) - bed_corners(1,2)
      fls_corners(2,2) = (rhoi_rhow*H_corners(2,2)) - bed_corners(2,2)
    endif

    ! Exact early-out: the sub-IP flotation deficit r*max(h_ip, min_h_shelf) - bed_ip
    ! is built from bilinear interpolants of the corner values gathered above, and a
    ! bilinear attains its extrema at the cell corners. If the unclamped corner
    ! deficits r*h - bed are all positive, the deficit is positive at every sub-IP
    ! (the min_h_shelf clamp can only raise it) and the sampled fraction is exactly
    ! 1; if they are all non-positive AND the clamped branch r*min_h_shelf - bed is
    ! non-positive at every corner, the deficit is non-positive at every sub-IP and
    ! the fraction is exactly 0. Both results are bitwise identical to sampling, so
    ! the 4*N^2 sub-sampling is confined to cells whose corner deficits mix signs
    ! (the grounding-line band).
    d_min = min(min(fls_corners(1,1), fls_corners(2,1)), &
                min(fls_corners(1,2), fls_corners(2,2)))
    d_max = max(max(fls_corners(1,1), fls_corners(2,1)), &
                max(fls_corners(1,2), fls_corners(2,2)))
    bed_min = min(min(bed_corners(1,1), bed_corners(2,1)), &
                  min(bed_corners(1,2), bed_corners(2,2)))
    ! On the FV sub-element path the MIN_H_SHELF clamp is applied at cell centers before the
    ! interpolation, so there is no post-interpolation clamp that could raise a sub-point deficit and
    ! the plain sign test on the corner deficits is exact on its own.
    if (fv_sub) then
      if (d_min > 0.0) then
        CS%ground_frac(i,j) = 1.0
        CS%xi_basal(i,j,:,:) = 0.0
        cycle
      elseif (d_max <= 0.0) then
        CS%ground_frac(i,j) = 0.0
        CS%xi_basal(i,j,:,:) = 1.0
        cycle
      endif
    ! Exact early-outs (bilinear extrema at the corners), reproducing the strict grounding test.
    else
      if (d_min > 0.0) then
        CS%ground_frac(i,j) = 1.0
        CS%xi_basal(i,j,:,:) = 0.0
        cycle
      elseif ((d_max <= 0.0) .and. ((rhoi_rhow*CS%min_h_shelf) - bed_min <= 0.0)) then
        CS%ground_frac(i,j) = 0.0
        CS%xi_basal(i,j,:,:) = 1.0
        cycle
      endif
    endif

    if (CS%use_sep2) then
      ! SEP2: exact Jacobian-weighted grounded area fraction of the sub-element
      ! partition, so the diagnostic and gate agree with the friction geometry.
      fls_gf(1) = fls_corners(1,1) ; fls_gf(2) = fls_corners(2,1)
      fls_gf(3) = fls_corners(1,2) ; fls_gf(4) = fls_corners(2,2)
      call sep2_cell_qps(fls_gf, nqp_gf, beta_gf, wref_gf, qpg_gf)
      do tq=1,4
        do kq=1,nqp_gf(tq)
          b1_gf = beta_gf(1,kq,tq) ; b2_gf = beta_gf(2,kq,tq)
          b3_gf = beta_gf(3,kq,tq) ; b4_gf = beta_gf(4,kq,tq)
          mS_gf = b1_gf + b2_gf ; mN_gf = b3_gf + b4_gf
          mW_gf = b1_gf + b3_gf ; mE_gf = b2_gf + b4_gf
          a_gf = (G%dxCv(i,j-1) * mS_gf) + (G%dxCv(i,j) * mN_gf)
          d_gf = (G%dyCu(i-1,j) * mW_gf) + (G%dyCu(i,j) * mE_gf)
          vt_gf(kq) = wref_gf(kq,tq) * (a_gf * d_gf)
          vg_gf(kq) = merge(vt_gf(kq), 0.0, qpg_gf(kq,tq))
          ! Same quadrature weighted by each corner's basis function, giving the shape-function
          ! weighted grounded fraction seen by that corner rather than by the cell as a whole.
          do cq=1,4
            vxt_gf(cq,kq) = vt_gf(kq) * beta_gf(cq,kq,tq)
            vxg_gf(cq,kq) = merge(vxt_gf(cq,kq), 0.0, qpg_gf(kq,tq))
          enddo
        enddo
        ! Orbit-grouped QP sums: (2,3), (4,5) and (6,7) are reflection pairs.
        if (nqp_gf(tq) == 3) then
          pg_gf(tq) = vg_gf(1) + (vg_gf(2) + vg_gf(3))
          pt_gf(tq) = vt_gf(1) + (vt_gf(2) + vt_gf(3))
          do cq=1,4
            pxg_gf(cq,tq) = vxg_gf(cq,1) + (vxg_gf(cq,2) + vxg_gf(cq,3))
            pxt_gf(cq,tq) = vxt_gf(cq,1) + (vxt_gf(cq,2) + vxt_gf(cq,3))
          enddo
        else
          pg_gf(tq) = (vg_gf(1) + (vg_gf(2) + vg_gf(3))) + &
                      ((vg_gf(4) + vg_gf(5)) + (vg_gf(6) + vg_gf(7)))
          pt_gf(tq) = (vt_gf(1) + (vt_gf(2) + vt_gf(3))) + &
                      ((vt_gf(4) + vt_gf(5)) + (vt_gf(6) + vt_gf(7)))
          do cq=1,4
            pxg_gf(cq,tq) = (vxg_gf(cq,1) + (vxg_gf(cq,2) + vxg_gf(cq,3))) + &
                            ((vxg_gf(cq,4) + vxg_gf(cq,5)) + (vxg_gf(cq,6) + vxg_gf(cq,7)))
            pxt_gf(cq,tq) = (vxt_gf(cq,1) + (vxt_gf(cq,2) + vxt_gf(cq,3))) + &
                            ((vxt_gf(cq,4) + vxt_gf(cq,5)) + (vxt_gf(cq,6) + vxt_gf(cq,7)))
          enddo
        endif
      enddo
      ! Opposite-pair grouping is invariant under any rotation/reflection (S=1,E=2,N=3,W=4).
      w_ground = (pg_gf(1) + pg_gf(3)) + (pg_gf(2) + pg_gf(4))
      w_total  = (pt_gf(1) + pt_gf(3)) + (pt_gf(2) + pt_gf(4))
      CS%ground_frac(i,j) = w_ground / w_total
      ! Role-grouped cross-triangle reduction: each corner takes each of the roles (A, B, farA,
      ! farB) exactly once over the 4 triangles (S=1, E=2, N=3, W=4), so under a rotation the
      ! corner and the triangles move together and the operand order is preserved. Same
      ! construction as CG_action_sep2_basal, which is what makes the result rotation-invariant
      ! to the bit. Corner order is SW, SE, NW, NE, matching fls_gf.
      xg_gf(1) = (pxg_gf(1,1) + pxg_gf(1,4)) + (pxg_gf(1,2) + pxg_gf(1,3))
      xg_gf(2) = (pxg_gf(2,2) + pxg_gf(2,1)) + (pxg_gf(2,3) + pxg_gf(2,4))
      xg_gf(3) = (pxg_gf(3,4) + pxg_gf(3,3)) + (pxg_gf(3,1) + pxg_gf(3,2))
      xg_gf(4) = (pxg_gf(4,3) + pxg_gf(4,2)) + (pxg_gf(4,4) + pxg_gf(4,1))
      xt_gf(1) = (pxt_gf(1,1) + pxt_gf(1,4)) + (pxt_gf(1,2) + pxt_gf(1,3))
      xt_gf(2) = (pxt_gf(2,2) + pxt_gf(2,1)) + (pxt_gf(2,3) + pxt_gf(2,4))
      xt_gf(3) = (pxt_gf(3,4) + pxt_gf(3,3)) + (pxt_gf(3,1) + pxt_gf(3,2))
      xt_gf(4) = (pxt_gf(4,3) + pxt_gf(4,2)) + (pxt_gf(4,4) + pxt_gf(4,1))
      CS%xi_basal(i,j,1,1) = 1.0 - (xg_gf(1) / xt_gf(1))
      CS%xi_basal(i,j,2,1) = 1.0 - (xg_gf(2) / xt_gf(2))
      CS%xi_basal(i,j,1,2) = 1.0 - (xg_gf(3) / xt_gf(3))
      CS%xi_basal(i,j,2,2) = 1.0 - (xg_gf(4) / xt_gf(4))
      cycle
    endif

    n_grounded = 0
    xi_num(:,:) = 0.0 ; xi_den(:,:) = 0.0
    do jsub=1,CS%n_sub_regularize ; do isub=1,CS%n_sub_regularize
      do jq=1,2 ; do iq=1,2
        if (fv_sub) then
          ! Interpolate the deficit itself, with the same bilinear interpolant the SEP3 sub-point
          ! test uses, so the grounded fraction, the friction, and the surface kink all key off one
          ! field. No post-interpolation clamp: MIN_H_SHELF was applied at the cell centers.
          g_ip = ((CS%Phisub(iq,jq,isub,jsub,1,1)*fls_corners(1,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,2,2)*fls_corners(2,2))) + &
                 ((CS%Phisub(iq,jq,isub,jsub,2,1)*fls_corners(2,1)) + &
                  (CS%Phisub(iq,jq,isub,jsub,1,2)*fls_corners(1,2)))
        else
          if (CS%use_DG_thickness) then
            h_ip = ((CS%Phisub(iq,jq,isub,jsub,1,1)*CS%h_nodal(i,j,1,1)) + &
                    (CS%Phisub(iq,jq,isub,jsub,2,2)*CS%h_nodal(i,j,2,2))) + &
                   ((CS%Phisub(iq,jq,isub,jsub,2,1)*CS%h_nodal(i,j,2,1)) + &
                    (CS%Phisub(iq,jq,isub,jsub,1,2)*CS%h_nodal(i,j,1,2)))
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

          g_ip = rhoi_rhow * h_ip - bed_ip
        endif
        if (g_ip > 0.0) n_grounded = n_grounded + 1
        ! Shape-function weighted floating fraction per corner: the same sub-point flotation test,
        ! accumulated against each corner's basis function instead of counted. Normalising by the
        ! same weights makes xi dimensionless, bounded by [0,1], and free of the cell metric.
        do ib=1,2 ; do ia=1,2
          wq_gf = CS%Phisub(iq,jq,isub,jsub,ia,ib)
          xi_den(ia,ib) = xi_den(ia,ib) + wq_gf
          if (g_ip <= 0.0) xi_num(ia,ib) = xi_num(ia,ib) + wq_gf
        enddo ; enddo
      enddo ; enddo
    enddo ; enddo

    CS%ground_frac(i,j) = real(n_grounded) / real(n_total)
    do ib=1,2 ; do ia=1,2
      CS%xi_basal(i,j,ia,ib) = xi_num(ia,ib) / xi_den(ia,ib)
    enddo ; enddo
  enddo ; enddo

  call pass_var(CS%ground_frac, G%Domain, complete=.true.)
  ! xi is read at every corner of every cell in the data domain by the averaged and sub-element
  ! nodal source operators, so its halo has to be current before the next advect step.
  call pass_corner_field(CS%xi_basal, G)

end subroutine compute_ground_frac

!> Grounded-area fraction of a rectangular region from the flotation function at its four
!! corners, via the analytic bilinear-flotation area integral of the quadrant grounding-line
!! parameterization. Ported from CISM's glissade_grounding_line.F90::compute_grounded_fraction
!! (W. Lipscomb, Los Alamos National Laboratory; CISM is LGPL), the scheme documented in
!! Leguy, Lipscomb & Asay-Davis (2021, The Cryosphere 15:3229-3253, sec. 2.2; orig. Leguy et al.
!! 2014). Corners are ordered counter-clockwise from the southwest. The flotation function is
!! positive where floating and non-positive where grounded, with the grounding line at the
!! contour f = 0; the routine returns the fraction of the region where f <= 0.
subroutine gl_quadrant_grounded_frac(f_flot, frac)
  real, dimension(4), intent(in) :: f_flot !< Flotation function at the 4 CCW corners (>0 floating) [Z ~> m]
  real,               intent(out) :: frac  !< Grounded area fraction of the region [nondim]

  real :: a, b, c, d           ! Coefficients of the bilinear f(x,y) = a + b*x + c*y + d*x*y on the unit square
  real :: f1, f2, f3, f4       ! Corner values after rotation to a canonical orientation
  real :: f_corner             ! Area of a single triangular/curved corner region
  real :: f_corner1, f_corner2 ! Areas of the two corner regions in the diagonal case
  real :: f_trapezoid          ! Floating-side trapezoid area in the adjacent case
  real :: var                  ! Sign selector for the diagonal case
  logical, dimension(4) :: cfloat, logvar
  logical :: adjacent          ! True if the two like corners are edge-adjacent (vs diagonal)
  logical :: rotated           ! True if the diagonal case was rotated 90 degrees before integrating
  integer :: nc, nfloat
  real, parameter :: eps06 = 1.0e-6 ! Guards the small-curvature (d->0) branches

  nfloat = 0
  do nc=1,4
    cfloat(nc) = (f_flot(nc) > 0.0)
    if (cfloat(nc)) nfloat = nfloat + 1
  enddo

  if (nfloat == 0) then
    frac = 1.0 ; return
  elseif (nfloat == 4 .or. minval(f_flot) == 0.0) then
    ! All floating, or grounded only where f == 0 exactly: treat as fully floating.
    frac = 0.0 ; return
  endif

  if (nfloat == 1 .or. nfloat == 3) then
    ! One corner is unlike the other three. Rotate it to the southwest (corner 1).
    if (nfloat == 1) then
      logvar(:) = cfloat(:)
    else
      logvar(:) = .not. cfloat(:)
    endif
    if (logvar(1)) then       ! no rotation
      f1 = f_flot(1) ; f2 = f_flot(2) ; f3 = f_flot(3) ; f4 = f_flot(4)
    elseif (logvar(2)) then   ! rotate 90 degrees
      f4 = f_flot(1) ; f1 = f_flot(2) ; f2 = f_flot(3) ; f3 = f_flot(4)
    elseif (logvar(3)) then   ! rotate 180 degrees
      f3 = f_flot(1) ; f4 = f_flot(2) ; f1 = f_flot(3) ; f2 = f_flot(4)
    else                      ! rotate 270 degrees
      f2 = f_flot(1) ; f3 = f_flot(2) ; f4 = f_flot(3) ; f1 = f_flot(4)
    endif
    a = f1 ; b = f2 - f1 ; c = f4 - f1 ; d = (f1 + f3) - (f2 + f4)
    ! Area of the corner region (floating if nfloat==1, grounded if nfloat==3):
    !   d /= 0: [(bc - ad) ln|1 - ad/(bc)| + ad] / d^2 ;  d -> 0: a^2 / (2 b c)
    if (abs((a*d)/(b*c)) > eps06) then
      f_corner = ((b*c - a*d) * log(abs(1.0 - (a*d)/(b*c))) + a*d) / (d*d)
    else
      f_corner = (a*a) / (2.0*b*c)
    endif
    if (nfloat == 1) then  ! f_corner is the floating area
      frac = 1.0 - f_corner
    else                   ! f_corner is the grounded area
      frac = f_corner
    endif

  else  ! nfloat == 2
    if (cfloat(1) .and. cfloat(2)) then       ! two floating corners adjacent; no rotation
      adjacent = .true. ; f1 = f_flot(1) ; f2 = f_flot(2) ; f3 = f_flot(3) ; f4 = f_flot(4)
    elseif (cfloat(2) .and. cfloat(3)) then   ! rotate 90 degrees
      adjacent = .true. ; f4 = f_flot(1) ; f1 = f_flot(2) ; f2 = f_flot(3) ; f3 = f_flot(4)
    elseif (cfloat(3) .and. cfloat(4)) then   ! rotate 180 degrees
      adjacent = .true. ; f3 = f_flot(1) ; f4 = f_flot(2) ; f1 = f_flot(3) ; f2 = f_flot(4)
    elseif (cfloat(4) .and. cfloat(1)) then   ! rotate 270 degrees
      adjacent = .true. ; f2 = f_flot(1) ; f3 = f_flot(2) ; f4 = f_flot(3) ; f1 = f_flot(4)
    else                                      ! two floating corners diagonally opposite
      adjacent = .false.
      var = f_flot(2)*f_flot(4) - f_flot(1)*f_flot(3)
      if (var >= 0.0) then
        f1 = f_flot(1) ; f2 = f_flot(2) ; f3 = f_flot(3) ; f4 = f_flot(4) ; rotated = .false.
      else
        f4 = f_flot(1) ; f1 = f_flot(2) ; f2 = f_flot(3) ; f3 = f_flot(4) ; rotated = .true.
      endif
    endif
    a = f1 ; b = f2 - f1 ; c = f4 - f1 ; d = (f1 + f3) - (f2 + f4)
    if (adjacent) then
      ! Floating-side trapezoid area:
      !   d /= 0: [(bc - ad) ln(1 + d/c) - bd] / d^2 ;  d -> 0: -(2a + b) / (2c)
      if (abs(d/c) > eps06) then
        f_trapezoid = ((b*c - a*d) * log(1.0 + d/c) - b*d) / (d*d)
      else
        f_trapezoid = -(2.0*a + b) / (2.0*c)
      endif
      frac = 1.0 - f_trapezoid
    else
      ! Two opposite corner regions; lower-left integral plus upper-right integral.
      if (abs(b*c - a*d) > eps06) then
        f_corner1 = ((b*c - a*d) * log(abs(1.0 - (a*d)/(b*c))) + a*d) / (d*d)
        f_corner2 = ((b*c - a*d) * log(abs((b*c - a*d)/((b+d)*(c+d)))) + d*((a+b)+(c+d))) / (d*d)
      else
        f_corner1 = (a*a) / (b*c)
        f_corner2 = ((a+b)*(a+c)) / (b*c)
      endif
      if (f_flot(1) > 0.0) then  ! southwest corner floating
        if (rotated) then ; frac = f_corner1 + f_corner2
        else              ; frac = 1.0 - (f_corner1 + f_corner2) ; endif
      else                       ! southwest corner grounded
        if (rotated) then ; frac = 1.0 - (f_corner1 + f_corner2)
        else              ; frac = f_corner1 + f_corner2 ; endif
      endif
    endif
  endif

  ! Guard against small round-off excursions outside the physical range.
  frac = min(max(frac, 0.0), 1.0)

end subroutine gl_quadrant_grounded_frac

!> Compute the analytic grounded ice fraction at B-grid nodes (CS%f_ground_node, for basal
!! friction) and at cell centers (CS%f_ground_cell, for the driving-stress surface blend) using
!! the quadrant grounding-line parameterization of Leguy et al. (2021). Each cell is split into
!! four quadrants; the grounded area of every quadrant is integrated analytically by
!! gl_quadrant_grounded_frac and the shared quadrant areas are summed to the surrounding node
!! (averaged over its 4 quadrants) and to the host cell (averaged over its 4 in-cell quadrants),
!! so the node and cell grounded areas stay mutually consistent. Cell-mean thickness (h_shelf)
!! and bed (bed_elev) are used, so this runs identically under FV or DG thickness advection.
subroutine compute_gl_quadrant_fractions(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure describing the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure used by the ice shelf

  real, dimension(SZDI_(G),SZDJ_(G)) :: f_flot ! Cell-center flotation function (>0 floating) [Z ~> m]
  logical, dimension(SZDI_(G),SZDJ_(G)) :: ice_cell ! True where ice is present (grounded or floating)
  real, allocatable, dimension(:,:,:) :: fgq    ! Grounded fraction of the 4 quadrants around each node [nondim]
  real, dimension(4) :: fv                       ! Flotation at the 4 CCW corners of one quadrant [Z ~> m]
  real :: rhoi_rhow                              ! Ice/ocean density ratio [nondim]
  real :: f_flot_land_min                        ! Minimum depth below sea level assigned to a land cell
                                                 ! in the LINEARB flotation function [Z ~> m]
  real :: f_flot_marine_min                      ! Minimum magnitude of f_flot in a marine cell under
                                                 ! the LINEARB flotation function [Z ~> m]
  real :: h_cell                                 ! Ice thickness used in the flotation function [Z ~> m]
  logical :: filled                              ! True once an ice-free cell has an ice neighbor to copy from
  logical :: vmask                               ! True if the node has at least one ice-covered neighbor cell
  integer :: i, j, ii, jj, isd, ied, jsd, jed, isc, iec, jsc, jec
  integer :: i_in, j_in, i_off, j_off, gisc, gjsc, giec, gjec

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  ! CISM's f_flotation_land_topg_min and f_flotation_marine_min (glissade_grounding_line).
  f_flot_land_min   = 50.0 * G%US%m_to_Z
  f_flot_marine_min = 1.0e-4 * G%US%m_to_Z
  i_off = G%idg_offset ; j_off = G%jdg_offset
  gisc = 1 ; gjsc = 1 ; giec = G%domain%niglobal ; gjec = G%domain%njglobal

  ! Cell-center flotation function (Leguy 2021 "linear" form = ocean cavity thickness): f > 0 where
  ! floating, f <= 0 where grounded, with the grounding line at f = 0. Sign matches the existing
  ! cell-center flotation test rhoi_rhow*h - bed (grounded when positive). Note CS%bed_elev is the
  ! depth of the bed below sea level, so it is CISM's -(topg - eus) and a land cell has bed_elev < 0.
  do j=jsd,jed ; do i=isd,ied
    ice_cell(i,j) = (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3)
    ! CISM HO_FLOTATION_FUNCTION_LINEARB: evaluate the same expression in every cell, so an ice-free
    ! cell reports its own bed rather than an extrapolated neighbor value, and no extrapolation pass
    ! is needed. Cells whose bed stands above sea level are land and are made strongly grounded with
    ! a floor on the assumed freeboard, and |f| is floored in marine cells for robustness.
    if (CS%bed_elev(i,j) <= 0.0) then   ! land: the bed stands at or above sea level
      f_flot(i,j) = min(CS%bed_elev(i,j), -f_flot_land_min)
    else
      h_cell = 0.0
      if (ice_cell(i,j)) h_cell = max(ISS%h_shelf(i,j), CS%min_h_shelf)
      f_flot(i,j) = CS%bed_elev(i,j) - rhoi_rhow * h_cell
      if (abs(f_flot(i,j)) < f_flot_marine_min) &
        f_flot(i,j) = sign(f_flot_marine_min, f_flot(i,j))
    endif
  enddo ; enddo

  ! Under LINEARB the halo across a solid (non-reentrant) wall would otherwise be read as real ocean:
  ! it holds bed_elev = 0 there (no neighbor PE fills it), which the expression above turns into a
  ! barely floating marine cell and which would then unground the wall nodes. Fill those cells by
  ! zero-gradient extension of the nearest in-domain cell, so a wall bounding grounded ice stays
  ! grounded and one bounding a shelf stays floating. The LINEAR branch does not need this: its
  ! ice-free cells carry no bed information at all and are handled by the extrapolation below.
  do j=jsd,jed ; do i=isd,ied
    i_in = i ; j_in = j
    if (.not. CS%reentrant_x) i_in = min(max(i+i_off, gisc), giec) - i_off
    if (.not. CS%reentrant_y) j_in = min(max(j+j_off, gjsc), gjec) - j_off
    if ((i_in /= i) .or. (j_in /= j)) then
      if ((i_in >= isd) .and. (i_in <= ied) .and. (j_in >= jsd) .and. (j_in <= jed)) &
        f_flot(i,j) = f_flot(i_in,j_in)
    endif
  enddo ; enddo

  call pass_var(f_flot, G%Domain)

  allocate(fgq(4,isd:ied,jsd:jed), source=0.0)
  CS%f_ground_node(:,:) = 0.0
  CS%f_ground_cell(:,:) = 0.0

  ! Per-node quadrant fractions. Node (I,J) is the NE corner of cell (i,j); the four cells around
  ! it are (i,j), (i+1,j), (i+1,j+1), (i,j+1) (the CISM vertex convention). Each quadrant's corner
  ! values are the (extrapolated) cell-center field interpolated to {cell center, two edge midpoints,
  ! node}. Only nodes with at least one ice-covered neighbor (vmask) are computed; nodes surrounded
  ! entirely by ice-free ocean stay floating (f_ground_node = 0), matching CISM's vmask gate.
  do j=jsd,jed-1 ; do i=isd,ied-1
    vmask = (ice_cell(i,j) .or. ice_cell(i+1,j)) .or. (ice_cell(i,j+1) .or. ice_cell(i+1,j+1))
    if (.not. vmask) cycle
    ! Quadrant 1: NE quarter of cell (i,j) (southwest of the node)
    fv(1) =        f_flot(i,j)
    fv(2) = 0.5 * (f_flot(i,j)   + f_flot(i+1,j))
    fv(3) = 0.25*((f_flot(i,j) + f_flot(i+1,j+1)) + (f_flot(i+1,j) + f_flot(i,j+1)))
    fv(4) = 0.5 * (f_flot(i,j)   + f_flot(i,j+1))
    call gl_quadrant_grounded_frac(fv, fgq(1,i,j))
    ! Quadrant 2: NW quarter of cell (i+1,j) (southeast of the node)
    fv(1) = 0.5 * (f_flot(i+1,j) + f_flot(i,j))
    fv(2) =        f_flot(i+1,j)
    fv(3) = 0.5 * (f_flot(i+1,j) + f_flot(i+1,j+1))
    fv(4) = 0.25*((f_flot(i,j) + f_flot(i+1,j+1)) + (f_flot(i+1,j) + f_flot(i,j+1)))
    call gl_quadrant_grounded_frac(fv, fgq(2,i,j))
    ! Quadrant 3: SW quarter of cell (i+1,j+1) (northeast of the node)
    fv(1) = 0.25*((f_flot(i,j) + f_flot(i+1,j+1)) + (f_flot(i+1,j) + f_flot(i,j+1)))
    fv(2) = 0.5 * (f_flot(i+1,j+1) + f_flot(i+1,j))
    fv(3) =        f_flot(i+1,j+1)
    fv(4) = 0.5 * (f_flot(i+1,j+1) + f_flot(i,j+1))
    call gl_quadrant_grounded_frac(fv, fgq(3,i,j))
    ! Quadrant 4: SE quarter of cell (i,j+1) (northwest of the node)
    fv(1) = 0.5 * (f_flot(i,j+1) + f_flot(i,j))
    fv(2) = 0.25*((f_flot(i,j) + f_flot(i+1,j+1)) + (f_flot(i+1,j) + f_flot(i,j+1)))
    fv(3) = 0.5 * (f_flot(i,j+1) + f_flot(i+1,j+1))
    fv(4) =        f_flot(i,j+1)
    call gl_quadrant_grounded_frac(fv, fgq(4,i,j))

    ! Quadrants 1..4 run SW, SE, NE, NW around the node, so a quarter turn permutes them
    ! cyclically. Pairing opposite quadrants makes the sum bit-for-bit rotation-invariant:
    ! the shift maps the two partial sums onto each other, and their addition commutes.
    ! Pairing adjacent quadrants would instead regroup the four addends and change the
    ! rounding. Same convention as the SEP2 and DG reductions elsewhere in this module.
    CS%f_ground_node(i,j) = 0.25*((fgq(1,i,j) + fgq(3,i,j)) + (fgq(2,i,j) + fgq(4,i,j)))
  enddo ; enddo

  ! Per-cell grounded fraction = mean of the four in-cell quadrants, taken from the node sums
  ! above so the cell and node grounded areas are built from the same quadrant areas. Cell (i,j)
  ! collects quadrant 3 of node (i-1,j-1), quadrant 4 of node (i,j-1), quadrant 1 of node (i,j),
  ! and quadrant 2 of node (i-1,j).
  do j=jsd+1,jed-1 ; do i=isd+1,ied-1
    CS%f_ground_cell(i,j) = 0.25*((fgq(3,i-1,j-1) + fgq(1,i,j)) + (fgq(4,i,j-1) + fgq(2,i-1,j)))
    ! Each quadrant is the quarter of the cell adjacent to one corner, so the same four numbers
    ! that average to the cell grounded fraction are, individually, the nodal grounded fractions.
    ! The floating fraction is their complement. No extra integration is needed.
    CS%xi_basal(i,j,1,1) = 1.0 - fgq(3,i-1,j-1)
    CS%xi_basal(i,j,2,1) = 1.0 - fgq(4,i,  j-1)
    CS%xi_basal(i,j,1,2) = 1.0 - fgq(2,i-1,j  )
    CS%xi_basal(i,j,2,2) = 1.0 - fgq(1,i,  j  )
  enddo ; enddo

  deallocate(fgq)
  call pass_corner_field(CS%xi_basal, G)
  call pass_var(CS%f_ground_cell, G%Domain)
  call pass_var(CS%f_ground_node, G%Domain, position=CORNER)

end subroutine compute_gl_quadrant_fractions


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

!> Finite-volume driving stress integrated over the sub-element grounding-line partition
!! (I2N_TAUD). Every quadrature point lies strictly on one side of the sub-element
!! grounding line and takes that side's surface-slope branch, so no point straddles the kink.
!!
!! The surface is reconstructed from the same two corner fields the friction uses:
!!
!!   S      = (1-r)*H + max(fls, 0)
!!   grad S = (1-r)*grad H  +  { grad fls   grounded
!!                             { 0          floating
!!
!! which reproduces grad H - grad bed where grounded and (1-r)*grad H where floating, while placing
!! the kink exactly on the fls = 0 contour that the partition cut on. The bed never appears.
!!
!! The thickness multiplying the slope is the cell mean, not an interpolant, so rho*g*H is constant
!! over the cell. That is also what makes the strong form used here identical to a surface-form
!! integration by parts: with H constant, rho*g*H*grad(S) = grad(rho*g*H*S) and the integrand
!! rho*g*H*S is continuous across the grounding line (S_grounded - S_floating = fls = 0 there), so
!! that IBP carries no internal grounding-line contour integral. The *pressure*-form IBP used by the
!! DG path does, because 1/2*rho*g*h^2 step-jumps by r/2*rho*g*h^2 at flotation; that is why the
!! calving-front Neumann term below, which is exactly that pressure, is kept as a separate face term
!! rather than folded into the volume integral.
!!
!! Under SEP2 the kinked term uses the P1-on-the-fan gradient, constant on each parent triangle,
!! because the SEP2 cut is the zero contour of that same P1 interpolant; using the bilinear gradient
!! there would put the kink on a slightly different curve and make S jump across the cut. Under SEP3
!! the sub-point flotation test is bilinear, so the bilinear gradient is the consistent choice. The
!! smooth term (1-r)*grad H is bilinear in both cases -- it carries no kink.
subroutine calc_shelf_driving_stress_i2n(CS, ISS, G, US, taudx, taudy)
  type(ice_shelf_dyn_CS), intent(in)   :: CS  !< A pointer to the ice shelf control structure
  type(ice_shelf_state), intent(in)    :: ISS !< A structure describing the ice-shelf state
  type(ocean_grid_type), intent(inout) :: G   !< The grid structure used by the ice shelf.
  type(unit_scale_type), intent(in)    :: US  !< A structure containing unit conversion factors
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudx !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: taudy !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]

  real, dimension(SZDIB_(G),SZDJB_(G),4) :: taudx_b, taudy_b ! Per-element corner contributions,
                       ! diagonal-pair summed for rotation invariance [R L3 Z T-2 ~> kg m s-2]
  real, dimension(4)     :: hc     ! Corner thickness, flattened SW,SE,NW,NE [Z ~> m]
  real, dimension(4)     :: fls    ! Corner flotation deficit, flattened SW,SE,NW,NE [Z ~> m]
  integer, dimension(4)  :: nqp    ! SEP2 QPs per parent triangle
  real, dimension(4,7,4) :: beta   ! SEP2 corner-basis weights per (corner, QP, triangle) [nondim]
  real, dimension(7,4)   :: wref   ! SEP2 reference measure per (QP, triangle) [nondim]
  logical, dimension(7,4) :: qpg   ! SEP2 grounded state per (QP, triangle)
  real, dimension(4,7)   :: valx, valy ! Per-QP nodal contributions [R L3 Z T-2 ~> kg m s-2]
  real, dimension(4,4)   :: px, py ! Per-(corner, triangle) partial sums [R L3 Z T-2 ~> kg m s-2]
  real, dimension(4)     :: dfxi, dfeta ! Per-triangle P1 gradient of fls in reference space [Z ~> m]
  real, dimension(4)     :: dhxi, dheta ! Per-triangle P1 gradient of H in reference space [Z ~> m]
  real, dimension(7) :: vsx, vsy ! Per-QP weight-scaled surface slopes [Z L ~> m2 m-1]
  real, dimension(7) :: vw       ! Per-QP quadrature weights [L2 ~> m2]
  real, dimension(4) :: psx, psy ! Per-triangle weighted-slope sums [Z L ~> m2 m-1]
  real, dimension(4) :: pw       ! Per-triangle weight sums [L2 ~> m2]
  real :: w_total           ! Total quadrature weight over the cell [L2 ~> m2]
  logical :: calc_slope_diag ! True if the surface-slope diagnostics are registered
  real :: b1, b2, b3, b4    ! Corner-basis weights at the QP [nondim]
  real :: mS, mN, mW, mE    ! Marginal sums: interpolation weights of the 4 cell edges [nondim]
  real :: a, d              ! Interpolated cell-edge spacings at the QP [L ~> m]
  real :: weight            ! Quadrature weight [L2 ~> m2]
  real :: dhdx_gp, dhdy_gp  ! Corner-field thickness gradients at the QP [Z L-1 ~> nondim]
  real :: dfdx_gp, dfdy_gp  ! Flotation-deficit gradients at the QP [Z L-1 ~> nondim]
  real :: dsdx_gp, dsdy_gp  ! Surface gradients at the QP [Z L-1 ~> nondim]
  real :: fx_gp, fy_gp      ! Driving-stress integrand at the QP [R L Z T-2 ~> kg m-1 s-2]
  real :: rho, rhow, rhoi_rhow ! Ice and ocean densities [R ~> kg m-3] and their ratio [nondim]
  real :: grav              ! The gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: He                ! Cell-mean ice thickness in the driving-stress prefactor [Z ~> m]
  real :: rgHe              ! rho*g*He, constant over the cell [R L2 Z T-2 ~> kg m-1 s-2]
  real :: smag, scale       ! Slope magnitude and MAX_SURFACE_SLOPE scaling [nondim]
  real :: neumann_val       ! Lateral-pressure boundary term [R Z L2 T-2 ~> kg s-2]
  real :: dxS, dxN, dyW, dyE ! Cell edge spacings [L ~> m]
  integer :: i, j, isc, iec, jsc, jec, t, k, c
  integer :: i_off, j_off, gisc, gjsc, giec, gjec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  i_off = G%idg_offset ; j_off = G%jdg_offset
  gisc = 1 ; gjsc = 1 ; giec = G%domain%niglobal ; gjec = G%domain%njglobal
  rho = CS%density_ice ; rhow = CS%density_ocean_avg ; grav = CS%g_Earth
  rhoi_rhow = rho/rhow

  taudx_b(:,:,:) = 0.0 ; taudy_b(:,:,:) = 0.0

  ! Surface-slope diagnostics: the quadrature-weighted cell mean of the same per-QP slope that is
  ! integrated above, so the diagnostic reports the branch-resolved slope the driving stress used
  ! rather than a separate reconstruction. Cells with no ice, or with an incomplete corner stencil,
  ! keep zero.
  calc_slope_diag = (CS%id_sx_shelf > 0 .or. CS%id_sy_shelf > 0 .or. CS%id_surf_slope_mag_shelf > 0)
  if (calc_slope_diag) then
    do j=jsc,jec ; do i=isc,iec
      CS%sx_shelf(i,j) = 0.0 ; CS%sy_shelf(i,j) = 0.0
    enddo ; enddo
  endif

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if (ISS%hmask(i,j) /= 1 .and. ISS%hmask(i,j) /= 3) cycle
    if (.not. (CS%corner_valid(I-1,J-1) .and. CS%corner_valid(I,J-1) .and. &
               CS%corner_valid(I-1,J  ) .and. CS%corner_valid(I,J  ))) cycle

    hc(1) = CS%H_corner(I-1,J-1) ; hc(2) = CS%H_corner(I,J-1)
    hc(3) = CS%H_corner(I-1,J  ) ; hc(4) = CS%H_corner(I,J  )
    fls(1) = CS%fls_corner(I-1,J-1) ; fls(2) = CS%fls_corner(I,J-1)
    fls(3) = CS%fls_corner(I-1,J  ) ; fls(4) = CS%fls_corner(I,J  )

    He = max(ISS%h_shelf(i,j), CS%min_h_shelf)
    rgHe = (rho*grav) * He
    dxS = G%dxCv(i,j-1) ; dxN = G%dxCv(i,j)
    dyW = G%dyCu(i-1,j) ; dyE = G%dyCu(i,j)

    call sep2_cell_qps(fls, nqp, beta, wref, qpg)

    ! Per-triangle P1 gradients of BOTH corner fields. Both terms of S must use the same operator
    ! or the grounded branch does not telescope to grad(h) - grad(bed); see sep2_fan_gradient.
    call sep2_fan_gradient(fls, dfxi, dfeta)
    call sep2_fan_gradient(hc,  dhxi, dheta)

    do t=1,4
      do k=1,nqp(t)
        b1 = beta(1,k,t) ; b2 = beta(2,k,t) ; b3 = beta(3,k,t) ; b4 = beta(4,k,t)
        mS = b1 + b2 ; mN = b3 + b4 ; mW = b1 + b3 ; mE = b2 + b4
        a = (dxS * mS) + (dxN * mN)
        d = (dyW * mW) + (dyE * mE)
        weight = wref(k,t) * (a * d)

        ! S = (1-r)*H + max(fls,0). Smooth part, no kink; P1 to match the kinked part below.
        dhdx_gp = dhxi(t) / a ; dhdy_gp = dheta(t) / d
        ! Kinked part: zero on the floating side, so grad S telescopes to P[h] - P[bed] where
        ! grounded and (1-r)*P[h] where floating, with the kink exactly on the cut.
        if (qpg(k,t)) then
          dfdx_gp = dfxi(t) / a ; dfdy_gp = dfeta(t) / d
        else
          dfdx_gp = 0.0 ; dfdy_gp = 0.0
        endif
        dsdx_gp = ((1.0 - rhoi_rhow) * dhdx_gp) + dfdx_gp
        dsdy_gp = ((1.0 - rhoi_rhow) * dhdy_gp) + dfdy_gp

        if (CS%max_surface_slope > 0) then
          smag = sqrt((dsdx_gp**2) + (dsdy_gp**2))
          scale = CS%max_surface_slope / max(smag, CS%max_surface_slope)
          dsdx_gp = scale*dsdx_gp ; dsdy_gp = scale*dsdy_gp
        endif

        fx_gp = -rgHe * dsdx_gp
        fy_gp = -rgHe * dsdy_gp
        do c=1,4
          valx(c,k) = (weight * beta(c,k,t)) * fx_gp
          valy(c,k) = (weight * beta(c,k,t)) * fy_gp
        enddo
        if (calc_slope_diag) then
          vw(k) = weight
          vsx(k) = dsdx_gp * weight
          vsy(k) = dsdy_gp * weight
        endif
      enddo

      ! Orbit-grouped QP sums: (2,3), (4,5) and (6,7) are reflection pairs.
      if (nqp(t) == 3) then
        do c=1,4
          px(c,t) = valx(c,1) + (valx(c,2) + valx(c,3))
          py(c,t) = valy(c,1) + (valy(c,2) + valy(c,3))
        enddo
        if (calc_slope_diag) then
          psx(t) = vsx(1) + (vsx(2) + vsx(3))
          psy(t) = vsy(1) + (vsy(2) + vsy(3))
          pw(t)  = vw(1) + (vw(2) + vw(3))
        endif
      else
        do c=1,4
          px(c,t) = (valx(c,1) + (valx(c,2) + valx(c,3))) + &
                    ((valx(c,4) + valx(c,5)) + (valx(c,6) + valx(c,7)))
          py(c,t) = (valy(c,1) + (valy(c,2) + valy(c,3))) + &
                    ((valy(c,4) + valy(c,5)) + (valy(c,6) + valy(c,7)))
        enddo
        if (calc_slope_diag) then
          psx(t) = (vsx(1) + (vsx(2) + vsx(3))) + ((vsx(4) + vsx(5)) + (vsx(6) + vsx(7)))
          psy(t) = (vsy(1) + (vsy(2) + vsy(3))) + ((vsy(4) + vsy(5)) + (vsy(6) + vsy(7)))
          pw(t)  = (vw(1) + (vw(2) + vw(3))) + ((vw(4) + vw(5)) + (vw(6) + vw(7)))
        endif
      endif
    enddo

    if (calc_slope_diag) then
      if ((i >= isc) .and. (i <= iec) .and. (j >= jsc) .and. (j <= jec)) then
        ! Opposite-pair grouping (S,N) and (E,W) is invariant under the triangle permutations of
        ! any rotation or reflection (S=1, E=2, N=3, W=4).
        w_total = (pw(1) + pw(3)) + (pw(2) + pw(4))
        if (w_total > 0.0) then
          CS%sx_shelf(i,j) = ((psx(1) + psx(3)) + (psx(2) + psx(4))) / w_total
          CS%sy_shelf(i,j) = ((psy(1) + psy(3)) + (psy(2) + psy(4))) / w_total
        endif
      endif
    endif

    ! Role-grouped cross-triangle reduction (see CG_action_sep2_basal).
    taudx_b(I-1,J-1,4) = (px(1,1) + px(1,4)) + (px(1,2) + px(1,3))
    taudx_b(I  ,J-1,3) = (px(2,2) + px(2,1)) + (px(2,3) + px(2,4))
    taudx_b(I-1,J  ,2) = (px(3,4) + px(3,3)) + (px(3,1) + px(3,2))
    taudx_b(I  ,J  ,1) = (px(4,3) + px(4,2)) + (px(4,4) + px(4,1))
    taudy_b(I-1,J-1,4) = (py(1,1) + py(1,4)) + (py(1,2) + py(1,3))
    taudy_b(I  ,J-1,3) = (py(2,2) + py(2,1)) + (py(2,3) + py(2,4))
    taudy_b(I-1,J  ,2) = (py(3,4) + py(3,3)) + (py(3,1) + py(3,2))
    taudy_b(I  ,J  ,1) = (py(4,3) + py(4,2)) + (py(4,4) + py(4,1))
  enddo ; enddo

  do J=jsc-1,jec ; do I=isc-1,iec
    taudx(I,J) = taudx(I,J) + ((taudx_b(I,J,1)+taudx_b(I,J,4)) + (taudx_b(I,J,2)+taudx_b(I,J,3)))
    taudy(I,J) = taudy(I,J) + ((taudy_b(I,J,1)+taudy_b(I,J,4)) + (taudy_b(I,J,2)+taudy_b(I,J,3)))
  enddo ; enddo

  ! Lateral-pressure (Neumann) boundary conditions at calving fronts and stress faces, identical to
  ! calc_shelf_driving_stress. This is the pressure-form face term and is deliberately kept separate
  ! from the volume integral above (see the header note on the two integrations by parts).
  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then
      if (CS%ground_frac(i,j) == 1) then
        neumann_val = ((.5 * grav) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2 - &
                                      rhow * max(0.0, CS%bed_elev(i,j))**2))
      else
        neumann_val = (.5 * grav) * ((1-rho/rhow) * (rho * max(ISS%h_shelf(i,j),CS%min_h_shelf)**2))
      endif
      if ((CS%u_face_mask_bdry(I-1,j) == 2) .OR. &
        ((ISS%hmask(i-1,j) == 0 .OR. ISS%hmask(i-1,j) == 2) .AND. (CS%reentrant_x .OR. (i+i_off /= gisc)))) then
        taudx(I-1,J-1) = taudx(I-1,J-1) - .5 * G%dyCu(I-1,j) * neumann_val
        taudx(I-1,J) = taudx(I-1,J) - .5 * G%dyCu(I-1,j) * neumann_val
      endif
      if ((CS%u_face_mask_bdry(I,j) == 2) .OR. &
        ((ISS%hmask(i+1,j) == 0 .OR. ISS%hmask(i+1,j) == 2) .and. (CS%reentrant_x .OR. (i+i_off /= giec)))) then
        taudx(I,J-1) = taudx(I,J-1) + .5 * G%dyCu(I,j) * neumann_val
        taudx(I,J) = taudx(I,J) + .5 * G%dyCu(I,j) * neumann_val
      endif
      if ((CS%v_face_mask_bdry(i,J-1) == 2) .OR. &
        ((ISS%hmask(i,j-1) == 0 .OR. ISS%hmask(i,j-1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjsc)))) then
        taudy(I-1,J-1) = taudy(I-1,J-1) - .5 * G%dxCv(i,J-1) * neumann_val
        taudy(I,J-1) = taudy(I,J-1) - .5 * G%dxCv(i,J-1) * neumann_val
      endif
      if ((CS%v_face_mask_bdry(i,J) == 2) .OR. &
        ((ISS%hmask(i,j+1) == 0 .OR. ISS%hmask(i,j+1) == 2) .and. (CS%reentrant_y .OR. (j+j_off /= gjec)))) then
        taudy(I-1,J) = taudy(I-1,J) + .5 * G%dxCv(i,J) * neumann_val
        taudy(I,J) = taudy(I,J) + .5 * G%dxCv(i,J) * neumann_val
      endif
    endif
  enddo ; enddo

end subroutine calc_shelf_driving_stress_i2n

!> Pre-compute the dual-cell Q1 (Lagrange) interpolation weights used to carry the cell-centered
!! thickness and flotation deficit to B-grid corners for the FV_SUBGRID_GL_* paths. Node (I,J) is
!! the NE corner of cell (i,j); its four cells are (i,j), (i+1,j), (i,j+1), (i+1,j+1), stored in
!! CS%corner_wt in the order SW, SE, NW, NE.
!!
!! Each cell is weighted by the *opposite* cell's spacing -- that is, by the proximity of its
!! centroid to the node -- which is the separable Q1 interpolant on the dual cell (the box whose
!! four corners are the surrounding cell centers) and reproduces a linear field exactly at the node.
!! An area-weighted mean does not: it weights toward the larger cell, whose centroid is farther from
!! the node, and for a linear field on cells of width h1 and h2 it is in error by (h2-h1)/2. That
!! also rules out the lumped-mass conservative projection (Huth et al. 2021, JAMES,
!! 10.1029/2020MS002277, eq. 29) reduced to one point per cell: a Q1 corner basis evaluated at its
!! own element centroid is 1/4 for every corner, so the basis factor cancels and only the area
!! weight is left. That form is correct for its purpose -- many scattered particles, where A_p is the
!! material a particle represents -- but here nothing downstream conserves the corner fields (the
!! total driving force is conserved by the Q1 partition of unity in the assembly and by the cell-mean
!! thickness in the prefactor), so the interpolatory weights are the right ones. Control-volume
!! averages elsewhere (CS%area_node, the nodal C, lumped_corner_mass) correctly stay area-weighted.
!!
!! On a uniform Cartesian grid every weight is exactly 0.25, so this reduces bitwise to a plain
!! four-cell mean and only differs where the grid is stretched.
subroutine build_corner_lagrange_weights(CS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf

  real :: dxW, dxE  ! Mean cell width of the west and east columns at a node [L ~> m]
  real :: dyS, dyN  ! Mean cell height of the south and north rows at a node [L ~> m]
  real :: wxW, wxE  ! Lagrange weights of the west and east columns [nondim]
  real :: wyS, wyN  ! Lagrange weights of the south and north rows [nondim]
  integer :: i, j, isd, ied, jsd, jed

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%corner_wt(:,:,:) = 0.25
  do j=jsd,jed-1 ; do i=isd,ied-1
    dxW = 0.5 * (G%dxT(i,  j) + G%dxT(i,  j+1))
    dxE = 0.5 * (G%dxT(i+1,j) + G%dxT(i+1,j+1))
    dyS = 0.5 * (G%dyT(i,j  ) + G%dyT(i+1,j  ))
    dyN = 0.5 * (G%dyT(i,j+1) + G%dyT(i+1,j+1))
    if ((dxW + dxE) <= 0.0 .or. (dyS + dyN) <= 0.0) cycle
    ! Weight each column/row by the opposite one's spacing (proximity of the centroid to the node).
    wxW = dxE / (dxW + dxE) ; wxE = dxW / (dxW + dxE)
    wyS = dyN / (dyS + dyN) ; wyN = dyS / (dyS + dyN)
    CS%corner_wt(1,I,J) = wxW * wyS  ! SW cell (i,j)
    CS%corner_wt(2,I,J) = wxE * wyS  ! SE cell (i+1,j)
    CS%corner_wt(3,I,J) = wxW * wyN  ! NW cell (i,j+1)
    CS%corner_wt(4,I,J) = wxE * wyN  ! NE cell (i+1,j+1)
  enddo ; enddo

end subroutine build_corner_lagrange_weights

!> Build the two corner fields that drive every sub-element grounding-line decision in the FV
!! (non-DG) path: the ice thickness CS%H_corner and the flotation deficit CS%fls_corner = r*h - bed.
!!
!! Both are carried from cell centers with the same weights (CS%corner_wt) over the same cell set, so
!! the identity fls = r*H - bed survives the interpolation. That is what lets the surface be written
!! as S = (1-r)*H + max(fls,0), whose slope kink is the zero contour of fls -- the same contour the
!! sub-element partition cuts on, and the same one the friction tests -- so the driving stress, the
!! basal friction, the effective pressure, and the grounded fraction cannot disagree about where the
!! grounding line is. The bed does not appear again after this routine; where it is needed (the
!! effective pressure cap) it is recovered as bed = r*H - fls.
!!
!! Cell inclusion follows the "option 3" ice-margin rule of Lipscomb et al. (2019), transplanted from
!! edges (edge_ok, in calc_shelf_driving_stress_vertex) to cell inclusion, because this interpolation
!! replaces the eq.-14 nodal gradient that used to host it:
!!   - ice-covered cells: always included, with h clamped at MIN_H_SHELF *before* the interpolation
!!     so that the clamped h and the fls built from it stay mutually consistent;
!!   - ice-free land lying below the ice (a real terrestrial margin): included with h = 0 exactly,
!!     for which S = max(-bed,0) is the bare-ground elevation and the effective pressure is zero, so
!!     no special case is needed anywhere downstream;
!!   - ice-free land standing above the ice (a nunatak): excluded. Including it would give the
!!     adjacent ice a surface sloping off the rock and a spurious driving stress; a nunatak is a
!!     lateral boundary that drags on the ice, not a source of driving stress;
!!   - ice-free ocean: excluded. Including it at any positive thickness drives the corner flotation
!!     deficit strongly negative, which floats a *grounded* marine terminus -- the lateral load at a
!!     calving front is already supplied by the Neumann face term;
!!   - cells outside a non-reentrant computational boundary: excluded. The halo across a solid wall
!!     has bed_elev = 0 with no neighbor PE to fill it, which is not a bed of zero but an undefined
!!     one. Excluding it from both the numerator and the weight sum drops it cleanly, unlike setting
!!     an ice-free flotation function to 0, which instead places the cell exactly on the flotation
!!     contour and reads as grounded.
!! Because ice-free cells are excluded rather than filled, no extrapolation pass into ice-free cells
!! is needed.
subroutine build_corner_flotation_fields(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS  !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure describing the ice-shelf state
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure used by the ice shelf

  integer, parameter :: KIND_SKIP = 0 !< Cell contributes to no corner
  integer, parameter :: KIND_ICE  = 1 !< Ice-covered cell
  integer, parameter :: KIND_LAND = 2 !< Ice-free land, included only where it lies below the ice
  integer, dimension(SZDI_(G),SZDJ_(G)) :: ckind ! Per-cell classification, one of KIND_*
  real, dimension(SZDI_(G),SZDJ_(G)) :: h_c   ! Cell thickness entering the interpolation [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: fls_c ! Cell flotation deficit r*h_c - bed [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: S_c   ! Cell surface elevation (1-r)*h_c + max(fls_c,0) [Z ~> m]
  real    :: rhoi_rhow    ! Ice/ocean density ratio r [nondim]
  real    :: S_ice_max    ! Largest surface elevation among the ice cells at a node [Z ~> m]
  real    :: wsum         ! Sum of the weights of the included cells [nondim]
  real    :: whs(4), wfs(4), wws(4) ! Weighted thickness, deficit and weight of the 4 cells,
                          ! zero where the cell is excluded [Z ~> m], [Z ~> m], [nondim]
  integer :: kc(4)        ! Classification of the 4 cells around a node
  real    :: sc4(4)       ! Surface elevation of the 4 cells around a node [Z ~> m]
  logical :: have_ice     ! True if at least one of the 4 cells is ice-covered
  integer :: i, j, n, isd, ied, jsd, jed
  integer :: ii(4), jj(4) ! Tracer indices of the 4 cells around a node, ordered SW, SE, NW, NE
  integer :: i_off, j_off, gisc, gjsc, giec, gjec

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  i_off = G%idg_offset ; j_off = G%jdg_offset
  gisc = 1 ; gjsc = 1 ; giec = G%domain%niglobal ; gjec = G%domain%njglobal
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg

  ! Cell-centered preparation. The MIN_H_SHELF clamp is applied here, before the interpolation, so
  ! that fls_c = r*h_c - bed holds for the clamped thickness and therefore survives to the corners.
  ckind(:,:) = KIND_SKIP ; h_c(:,:) = 0.0 ; fls_c(:,:) = 0.0 ; S_c(:,:) = 0.0
  do j=jsd,jed ; do i=isd,ied
    ! Cells outside a non-reentrant computational boundary have no meaningful bed and are skipped.
    if (.not. CS%reentrant_x) then
      if ((i+i_off < gisc) .or. (i+i_off > giec)) cycle
    endif
    if (.not. CS%reentrant_y) then
      if ((j+j_off < gjsc) .or. (j+j_off > gjec)) cycle
    endif
    if (ISS%hmask(i,j) == 1 .or. ISS%hmask(i,j) == 3) then
      ckind(i,j) = KIND_ICE
      h_c(i,j) = max(ISS%h_shelf(i,j), CS%min_h_shelf)
    elseif (CS%bed_elev(i,j) <= 0.0) then
      ckind(i,j) = KIND_LAND
      h_c(i,j) = 0.0
    else
      cycle  ! ice-free ocean
    endif
    fls_c(i,j) = (rhoi_rhow * h_c(i,j)) - CS%bed_elev(i,j)
    S_c(i,j) = ((1.0 - rhoi_rhow) * h_c(i,j)) + max(fls_c(i,j), 0.0)
  enddo ; enddo

  CS%H_corner(:,:) = 0.0 ; CS%fls_corner(:,:) = 0.0 ; CS%corner_valid(:,:) = .false.

  do j=jsd,jed-1 ; do i=isd,ied-1
    ii(1) = i   ; jj(1) = j     ! SW
    ii(2) = i+1 ; jj(2) = j     ! SE
    ii(3) = i   ; jj(3) = j+1   ! NW
    ii(4) = i+1 ; jj(4) = j+1   ! NE
    do n=1,4
      kc(n) = ckind(ii(n),jj(n)) ; sc4(n) = S_c(ii(n),jj(n))
    enddo

    ! The ice-free-land test is against the highest ice surface among the cells sharing this node:
    ! "this ground stands above the ice" is the conservative reading of the option-3 nunatak rule.
    have_ice = .false. ; S_ice_max = 0.0
    do n=1,4
      if (kc(n) == KIND_ICE) then
        if (have_ice) then ; S_ice_max = max(S_ice_max, sc4(n))
        else ; S_ice_max = sc4(n) ; have_ice = .true. ; endif
      endif
    enddo
    if (.not. have_ice) cycle  ! no ice touches this node; its value is never used

    do n=1,4
      if ((kc(n) == KIND_ICE) .or. ((kc(n) == KIND_LAND) .and. (sc4(n) < S_ice_max))) then
        wws(n) = CS%corner_wt(n,I,J)
        whs(n) = wws(n) * h_c(ii(n),jj(n))
        wfs(n) = wws(n) * fls_c(ii(n),jj(n))
      else
        wws(n) = 0.0 ; whs(n) = 0.0 ; wfs(n) = 0.0
      endif
    enddo

    ! Diagonal-pair sums (SW+NE)+(SE+NW), so the corner fields are bitwise invariant under a
    ! 90-degree grid rotation, matching CG_action and lumped_corner_mass.
    wsum = (wws(1) + wws(4)) + (wws(2) + wws(3))
    if (wsum <= 0.0) cycle
    CS%H_corner(I,J)   = ((whs(1) + whs(4)) + (whs(2) + whs(3))) / wsum
    CS%fls_corner(I,J) = ((wfs(1) + wfs(4)) + (wfs(2) + wfs(3))) / wsum
    CS%corner_valid(I,J) = .true.
  enddo ; enddo

  call pass_var(CS%H_corner, G%domain, position=CORNER, complete=.false.)
  call pass_var(CS%fls_corner, G%domain, position=CORNER, complete=.true.)

end subroutine build_corner_flotation_fields

!> Coulomb fB parameter at a quadrature point from the effective pressure directly, for the
!! FV_SUBGRID_GL_* paths where N is formed from the corner flotation field rather than from a
!! thickness and a bed. Identical in form to compute_fB_local once N is in hand.
pure real function compute_fB_from_N(N_eff, C_basal, alpha_coulomb, CF_Max, CF_MinN, &
    CF_PostPeak, n_basal_fric)
  real, intent(in) :: N_eff         !< Effective pressure at the quadrature point [R Z L T-2 ~> Pa]
  real, intent(in) :: C_basal       !< Basal friction coefficient for this cell [R L Z T-2 (s m-1)^n]
  real, intent(in) :: alpha_coulomb !< Coulomb prefactor [nondim]
  real, intent(in) :: CF_Max        !< Coulomb friction maximum coefficient [nondim]
  real, intent(in) :: CF_MinN       !< Minimum Coulomb effective pressure [R Z L T-2 ~> Pa]
  real, intent(in) :: CF_PostPeak   !< Coulomb post-peak exponent q [nondim]
  real, intent(in) :: n_basal_fric  !< Friction sliding exponent m [nondim]

  real :: fN  ! Floored effective pressure [R Z L T-2 ~> Pa]

  fN = max(N_eff, CF_MinN)
  if (fN > 0.0) then
    compute_fB_from_N = alpha_coulomb * (C_basal / (CF_Max * fN))**(CF_PostPeak / n_basal_fric)
  else
    ! Zero effective pressure: the Coulomb drag is exactly zero, which this factorization can only
    ! express as an infinite fB. Signal it instead. Unreachable when CF_MinN > 0.
    compute_fB_from_N = FB_NO_COULOMB_DRAG
  endif
end function compute_fB_from_N

!> Effective pressure at a quadrature point from the sub-element flotation field:
!! N = rho_ocean*g*min(fls, r*H). Where the bed is below sea level this is the usual
!! rho_ice*g*(H - H_f); where it is above, min selects r*H and N reduces to the pure overburden
!! rho_ice*g*H, with no spurious water column over dry land. This is the same cap CISM applies by
!! clamping f_pattyn to [0,1], written without a branch, and it needs no bed field: the bed is
!! implicit in fls. Taking fls from the same interpolant that decided the quadrature point's
!! flotation state guarantees N >= 0 at every grounded point.
pure real function subgrid_effective_pressure(fls_qp, h_qp, rhoi_rhow, rho_ocean_g_LtoZ)
  real, intent(in) :: fls_qp    !< Flotation deficit r*h - bed at the quadrature point [Z ~> m]
  real, intent(in) :: h_qp      !< Ice thickness at the quadrature point [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio r [nondim]
  real, intent(in) :: rho_ocean_g_LtoZ !< US%L_to_Z * density_ocean_avg * g_Earth [R L Z-1 T-2]

  subgrid_effective_pressure = rho_ocean_g_LtoZ * min(fls_qp, rhoi_rhow * h_qp)
end function subgrid_effective_pressure

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

!> DG(1) version of interpolate_H_to_B: average the co-located corners of the
!! ice-covered cells at each B-grid node.
subroutine interpolate_H_to_B_DG(G, h_nodal, hmask, H_node, min_h_shelf)
  type(ocean_grid_type), intent(in) :: G  !< The grid structure used by the ice shelf.
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

  ! Cell (ic,jc) = (i-1+k,j-1+l) touches node (i,j) at its corner (3-k,3-l).
  do j=jsc-1,jec
    do i=isc-1,iec
      num_h = 0
      h_arr(:,:) = 0.0
      do l=1,2 ; jc=j-1+l ; do k=1,2 ; ic=i-1+k
        if (hmask(ic,jc) == 1.0 .or. hmask(ic,jc) == 3.0) then
          h_arr(k,l) = max(h_nodal(ic,jc,3-k,3-l), min_h_shelf)
          num_h = num_h + 1
        endif
      enddo ; enddo
      if (num_h > 0) then
        H_node(i,j) = ((h_arr(1,1)+h_arr(2,2))+(h_arr(1,2)+h_arr(2,1))) / num_h
      endif
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
  if (associated(CS%coef_prefactor_node)) deallocate(CS%coef_prefactor_node)
  if (associated(CS%fB_node)) deallocate(CS%fB_node)
  if (associated(CS%area_node)) deallocate(CS%area_node)
  deallocate(CS%OD_rt, CS%OD_av)
  if (associated(CS%H_node)) deallocate(CS%H_node)
  deallocate(CS%t_bdry_val, CS%bed_elev)
  if (associated(CS%bed_node)) deallocate(CS%bed_node)
  if (associated(CS%h_nodal)) deallocate(CS%h_nodal)
  if (associated(CS%Minv_nodal)) deallocate(CS%Minv_nodal)
  if (associated(CS%cell_mean_w)) deallocate(CS%cell_mean_w)
  if (associated(CS%h_source_rate)) deallocate(CS%h_source_rate)
  if (associated(CS%h_source_rate_bmb)) deallocate(CS%h_source_rate_bmb)
  if (associated(CS%xi_basal)) deallocate(CS%xi_basal)
  if (associated(CS%h_source_rate_last)) deallocate(CS%h_source_rate_last)
  if (associated(CS%phi_x_FV)) deallocate(CS%phi_x_FV)
  if (associated(CS%phi_y_FV)) deallocate(CS%phi_y_FV)
  if (associated(CS%dg_art_visc_coef_u)) deallocate(CS%dg_art_visc_coef_u)
  if (associated(CS%dg_art_visc_coef_v)) deallocate(CS%dg_art_visc_coef_v)
  if (associated(CS%dg_art_visc_nu_u)) deallocate(CS%dg_art_visc_nu_u)
  if (associated(CS%dg_art_visc_nu_v)) deallocate(CS%dg_art_visc_nu_v)
  if (associated(CS%dg_art_visc_cell_scale)) deallocate(CS%dg_art_visc_cell_scale)
  if (associated(CS%dg_slow_idle_face_u)) deallocate(CS%dg_slow_idle_face_u)
  if (associated(CS%dg_slow_idle_face_v)) deallocate(CS%dg_slow_idle_face_v)
  deallocate(CS%ground_frac, CS%ground_frac_rt)
  if (associated(CS%dg_damp_ahat)) deallocate(CS%dg_damp_ahat)
  if (associated(CS%dg_damp_spread)) deallocate(CS%dg_damp_spread)
  if (associated(CS%dg_damp_tend)) deallocate(CS%dg_damp_tend)
  if (associated(CS%dg_damp_gate)) deallocate(CS%dg_damp_gate)
  if (associated(CS%dg_damp_want)) deallocate(CS%dg_damp_want)
  if (associated(CS%dg_damp_got)) deallocate(CS%dg_damp_got)
  if (associated(CS%dg_tilt_relax_rate)) deallocate(CS%dg_tilt_relax_rate)
  if (associated(CS%dg_tilt_relax_tend)) deallocate(CS%dg_tilt_relax_tend)
  if (associated(CS%dg_tilt_relax_alt_x)) deallocate(CS%dg_tilt_relax_alt_x)
  if (associated(CS%dg_tilt_relax_alt_y)) deallocate(CS%dg_tilt_relax_alt_y)
  if (associated(CS%dg_tilt_relax_alt_w)) deallocate(CS%dg_tilt_relax_alt_w)
  if (associated(CS%dg_tilt_relax_state)) deallocate(CS%dg_tilt_relax_state)
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



!> Add the ice-front Neumann term int phi*(P_ice - P_ocean)*n dS on one face of a cell
!! to its two endpoint corners, with P_ice = rho*g*h^2/2 and P_ocean = rhow*g*d^2/2.
subroutine add_Neumann_face_DG(face_length, face_sign, &
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
    h_loc = ((1.0 - t_face)*h_A) + (t_face*h_B)
    b_loc = ((1.0 - t_face)*b_A) + (t_face*b_B)
    P_ice = 0.5 * grav * rho * h_loc**2
    d_ocean = max(0.0, min(b_loc, rhoi_rhow * h_loc))
    P_ocean = 0.5 * grav * rhow * d_ocean**2
    P_face = P_ice - P_ocean
    phi_A = 1.0 - t_face ; phi_B = t_face
    face_A = face_A + (face_sign * 0.5 * face_length * phi_A * P_face)
    face_B = face_B + (face_sign * 0.5 * face_length * phi_B * P_face)
  enddo
end subroutine add_Neumann_face_DG

!> Add this cell's half of the interior-face term -1/2 int phi*rho*g*{h}*[s]*n dS, which the
!! cell volume integral misses because s jumps across the face. Each side takes its own
!! flotation branch for s.
subroutine add_taud_edge_correction_DG(face_length, face_sign, &
    h_loc_A, h_loc_B, h_ngh_A, h_ngh_B, b_corner_A, b_corner_B, &
    rho, rhoi_rhow, grav, min_h_shelf, &
    xquad, face_A, face_B)
  real, intent(in)    :: face_length      !< Face length [L ~> m]
  real, intent(in)    :: face_sign        !< +1 or -1, outward unit-normal component [nondim]
  real, intent(in)    :: h_loc_A          !< Local-cell h_nodal at face endpoint A [Z ~> m]
  real, intent(in)    :: h_loc_B          !< Local-cell h_nodal at face endpoint B [Z ~> m]
  real, intent(in)    :: h_ngh_A          !< Neighbour-cell h_nodal at face endpoint A [Z ~> m]
  real, intent(in)    :: h_ngh_B          !< Neighbour-cell h_nodal at face endpoint B [Z ~> m]
  real, intent(in)    :: b_corner_A       !< bed depth at face endpoint A [Z ~> m]
  real, intent(in)    :: b_corner_B       !< bed depth at face endpoint B [Z ~> m]
  real, intent(in)    :: rho              !< Ice density [R ~> kg m-3]
  real, intent(in)    :: rhoi_rhow        !< Ice to ocean density ratio [nondim]
  real, intent(in)    :: grav             !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real, intent(in)    :: min_h_shelf      !< Lower clamp on h [Z ~> m]
  real, dimension(2), intent(in) :: xquad !< 2-point Gauss-Legendre nodes on [0,1] [nondim]
  real, intent(inout) :: face_A           !< Accumulator for corner A [R L3 Z T-2 ~> kg m s-2]
  real, intent(inout) :: face_B           !< Accumulator for corner B [R L3 Z T-2 ~> kg m s-2]

  real :: hL_A, hL_B, hN_A, hN_B, b_A, b_B
  real :: t_face, h_loc, h_ngh, b_loc
  real :: s_loc, s_ngh, h_avg, jump_factor
  real :: phi_A, phi_B
  integer :: gp_face

  hL_A = max(h_loc_A, min_h_shelf) ; hL_B = max(h_loc_B, min_h_shelf)
  hN_A = max(h_ngh_A, min_h_shelf) ; hN_B = max(h_ngh_B, min_h_shelf)
  b_A = b_corner_A ; b_B = b_corner_B

  do gp_face = 1, 2
    t_face = xquad(gp_face)
    h_loc = ((1.0 - t_face)*hL_A) + (t_face*hL_B)
    h_ngh = ((1.0 - t_face)*hN_A) + (t_face*hN_B)
    b_loc = ((1.0 - t_face)*b_A) + (t_face*b_B)

    if (rhoi_rhow * h_loc - b_loc > 0.0) then
      s_loc = h_loc - b_loc                          ! grounded
    else
      s_loc = (1.0 - rhoi_rhow) * h_loc              ! floating
    endif
    if (rhoi_rhow * h_ngh - b_loc > 0.0) then
      s_ngh = h_ngh - b_loc                          ! grounded
    else
      s_ngh = (1.0 - rhoi_rhow) * h_ngh              ! floating
    endif

    h_avg = 0.5 * (h_loc + h_ngh)
    jump_factor = rho * grav * h_avg * (s_loc - s_ngh)

    phi_A = 1.0 - t_face ; phi_B = t_face
    !0.25 is the Gauss weight times the cell's half share of the inter-cell edge integral
    face_A = face_A + (face_sign * 0.25 * face_length * phi_A * jump_factor)
    face_B = face_B + (face_sign * 0.25 * face_length * phi_B * jump_factor)
  enddo
end subroutine add_taud_edge_correction_DG

!> DG(1) driving stress: the volume integral int phi*(-rho*g*h*grad(s)) dA (2x2 Gauss, or
!! SEP2/SEP3 sub-element quadrature in grounding-line cells), plus the interior-face jump
!! term (add_taud_edge_correction_DG) and the ice-front Neumann term (add_Neumann_face_DG).
subroutine calc_shelf_driving_stress_DG(CS, ISS, G, taudx, taudy)
  type(ice_shelf_dyn_CS), intent(inout) :: CS !< The ice shelf dynamics control structure
  type(ice_shelf_state),  intent(in)    :: ISS !< A structure with elements that describe the ice-shelf state
  type(ocean_grid_type),  intent(inout) :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: taudx  !< X-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                          intent(inout) :: taudy  !< Y-direction driving stress at q-points [R L3 Z T-2 ~> kg m s-2]

  real :: rho        ! Ice density [R ~> kg m-3]
  real :: rhow       ! Reference ocean density [R ~> kg m-3]
  real :: rhoi_rhow  ! Ice/ocean density ratio [nondim]
  real :: grav       ! Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real :: h_gp                  ! Ice thickness at a qp [Z ~> m]
  real :: dhdx_gp, dhdy_gp      ! Thickness gradients in physical coords [Z L-1 ~> nondim]
  real :: bed_gp                ! Bed depth at a qp [Z ~> m]
  real :: dbdx_gp, dbdy_gp      ! Bed-depth gradients in physical coords [Z L-1 ~> nondim]
  real :: dbdx_ref, dbdy_ref    ! Bed-depth gradients in reference coords [Z ~> m]
  real :: dsdx_gp, dsdy_gp      ! Surface gradients in physical coords [nondim]
  real :: a_qp, d_qp            ! Per-qp interpolated cell-edge spacings [L ~> m]
  real :: weight                ! Per-qp quadrature weight including Jacobian [L2 ~> m2]
  real :: phi_val               ! Bilinear nodal basis value at a qp [nondim]
  real :: bed_corners(2,2)      ! Bed depth at the 4 B-grid corners of an element [Z ~> m]
  real :: dxCv_S, dxCv_N        ! Cell-edge spacings on south and north faces [L ~> m]
  real :: dyCu_W, dyCu_E        ! Cell-edge spacings on west and east faces [L ~> m]
  real :: fx_gp, fy_gp          ! Driving-force density at a qp [R Z L T-2 ~> kg m-1 s-2]

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

  real, dimension(SZDIB_(G),SZDJB_(G),4) :: taudx_b, taudy_b ! Node driving stress from each
                                           ! surrounding cell, 1=SW, 2=SE, 3=NW, 4=NE [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2) :: xquad
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

  do j=jsc-1,jec+1 ; do i=isc-1,iec+1
    if (ISS%hmask(i,j) /= 1 .and. ISS%hmask(i,j) /= 3) cycle

    bed_corners(1,1) = CS%bed_node(I-1,J-1)
    bed_corners(2,1) = CS%bed_node(I,J-1)
    bed_corners(1,2) = CS%bed_node(I-1,J)
    bed_corners(2,2) = CS%bed_node(I,J)

    dxCv_S = G%dxCv(i,J-1) ; dxCv_N = G%dxCv(i,J)
    dyCu_W = G%dyCu(I-1,j) ; dyCu_E = G%dyCu(I,j)

    ! Volume integral: sub-element quadrature where the cell straddles flotation.
    use_subgrid_cell = CS%GL_regularize .and. &
                       (CS%ground_frac(i,j) > 0.0) .and. (CS%ground_frac(i,j) < 1.0)
    if (use_subgrid_cell) then
      if (CS%use_sep2) then
        call calc_shelf_driving_stress_DG_sep2(CS, CS%h_nodal(i,j,:,:), bed_corners, &
            dxCv_S, dxCv_N, dyCu_W, dyCu_E, rho, rhoi_rhow, grav, vol_dx, vol_dy, &
            CS%sx_shelf(i,j), CS%sy_shelf(i,j), calc_slope_diag)
      else
        call calc_shelf_driving_stress_DG_subgrid(CS, CS%Phisub, &
            CS%h_nodal(i,j,:,:), bed_corners, &
            dxCv_S, dxCv_N, dyCu_W, dyCu_E, &
            rho, rhoi_rhow, grav, vol_dx, vol_dy, &
            CS%sx_shelf(i,j), CS%sy_shelf(i,j), calc_slope_diag)
      endif
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
        dhdx_gp = ( (((-xquad(3-jq))*CS%h_nodal(i,j,1,1)) + (( xquad(jq))  *CS%h_nodal(i,j,2,2))) + &
                    ((( xquad(3-jq))*CS%h_nodal(i,j,2,1)) + ((-xquad(jq))  *CS%h_nodal(i,j,1,2))) ) / a_qp
        dhdy_gp = ( (((-xquad(3-iq))*CS%h_nodal(i,j,1,1)) + (( xquad(iq))  *CS%h_nodal(i,j,2,2))) + &
                    (((-xquad(iq))  *CS%h_nodal(i,j,2,1)) + (( xquad(3-iq))*CS%h_nodal(i,j,1,2))) ) / d_qp

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

        if (CS%GL_couple) then
          is_grounded = (CS%ground_frac(i,j) >= 1.0)
        else
          is_grounded = (rhoi_rhow * h_gp - bed_gp > 0.0)
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

    ! Face terms: Neumann at ice fronts, the jump term at faces to ice. Walls add nothing.
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
      call add_Neumann_face_DG(G%dyCu(I-1,j), -1.0, &
        CS%h_nodal(i,j,1,1), CS%h_nodal(i,j,1,2), bed_corners(1,1), bed_corners(1,2), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dx_W_A, face_dx_W_B)
    elseif (ISS%hmask(i-1,j) == 1.0) then
      call add_taud_edge_correction_DG(G%dyCu(I-1,j), -1.0, &
        CS%h_nodal(i,j,1,1),   CS%h_nodal(i,j,1,2), &
        CS%h_nodal(i-1,j,2,1), CS%h_nodal(i-1,j,2,2), &
        bed_corners(1,1), bed_corners(1,2), &
        rho, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dx_W_A, face_dx_W_B)
    endif

    ! East face (n_x = +1).
    is_ext_bdry = ((CS%u_face_mask_bdry(I,j) == 2) .or. &
                  ((ISS%hmask(i+1,j) == 0 .or. ISS%hmask(i+1,j) == 2) .and. &
                   (CS%reentrant_x .or. (i+i_off /= giec))))
    if (is_ext_bdry) then
      call add_Neumann_face_DG(G%dyCu(I,j), +1.0, &
        CS%h_nodal(i,j,2,1), CS%h_nodal(i,j,2,2), bed_corners(2,1), bed_corners(2,2), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dx_E_A, face_dx_E_B)
    elseif (ISS%hmask(i+1,j) == 1.0) then
      call add_taud_edge_correction_DG(G%dyCu(I,j), +1.0, &
        CS%h_nodal(i,j,2,1),   CS%h_nodal(i,j,2,2), &
        CS%h_nodal(i+1,j,1,1), CS%h_nodal(i+1,j,1,2), &
        bed_corners(2,1), bed_corners(2,2), &
        rho, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dx_E_A, face_dx_E_B)
    endif

    ! South face (n_y = -1).
    is_ext_bdry = ((CS%v_face_mask_bdry(i,J-1) == 2) .or. &
                  ((ISS%hmask(i,j-1) == 0 .or. ISS%hmask(i,j-1) == 2) .and. &
                   (CS%reentrant_y .or. (j+j_off /= gjsc))))
    if (is_ext_bdry) then
      call add_Neumann_face_DG(G%dxCv(i,J-1), -1.0, &
        CS%h_nodal(i,j,1,1), CS%h_nodal(i,j,2,1), bed_corners(1,1), bed_corners(2,1), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dy_S_A, face_dy_S_B)
    elseif (ISS%hmask(i,j-1) == 1.0) then
      call add_taud_edge_correction_DG(G%dxCv(i,J-1), -1.0, &
        CS%h_nodal(i,j,1,1),   CS%h_nodal(i,j,2,1), &
        CS%h_nodal(i,j-1,1,2), CS%h_nodal(i,j-1,2,2), &
        bed_corners(1,1), bed_corners(2,1), &
        rho, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dy_S_A, face_dy_S_B)
    endif

    ! North face (n_y = +1).
    is_ext_bdry = ((CS%v_face_mask_bdry(i,J) == 2) .or. &
                  ((ISS%hmask(i,j+1) == 0 .or. ISS%hmask(i,j+1) == 2) .and. &
                   (CS%reentrant_y .or. (j+j_off /= gjec))))
    if (is_ext_bdry) then
      call add_Neumann_face_DG(G%dxCv(i,J), +1.0, &
        CS%h_nodal(i,j,1,2), CS%h_nodal(i,j,2,2), bed_corners(1,2), bed_corners(2,2), &
        CS%h_bdry_val(i,j), loc_is_bc, rho, rhow, rhoi_rhow, grav, CS%min_h_shelf, &
        xquad, face_dy_N_A, face_dy_N_B)
    elseif (ISS%hmask(i,j+1) == 1.0) then
      call add_taud_edge_correction_DG(G%dxCv(i,J), +1.0, &
        CS%h_nodal(i,j,1,2),   CS%h_nodal(i,j,2,2), &
        CS%h_nodal(i,j+1,1,1), CS%h_nodal(i,j+1,2,1), &
        bed_corners(1,2), bed_corners(2,2), &
        rho, rhoi_rhow, grav, CS%min_h_shelf, xquad, &
        face_dy_N_A, face_dy_N_B)
    endif

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

end subroutine calc_shelf_driving_stress_DG

!> DG(1) driving-stress volume integral over the nsub x nsub SEP3 sub-cells, with a
!! flotation test at each sub-quadrature point. The caller adds the face terms.
subroutine calc_shelf_driving_stress_DG_subgrid(CS, Phisub, &
    h_nodal_cell, bed_corners, &
    dxCv_S, dxCv_N, dyCu_W, dyCu_E, &
    rho, rhoi_rhow, grav, vol_dx, vol_dy, sx_shelf, sy_shelf, calc_slope_diag)
  type(ice_shelf_dyn_CS), intent(in) :: CS    !< Ice shelf control structure
  real, dimension(:,:,:,:,:,:), intent(in) :: Phisub !< Sub-grid quadrature weights [nondim]
  real, dimension(2,2), intent(in) :: h_nodal_cell !< Q1 nodal thickness at the 4 corners [Z ~> m]
  real, dimension(2,2), intent(in) :: bed_corners !< Bed depth at the 4 cell corners [Z ~> m]
  real, intent(in) :: dxCv_S         !< Cell x-length on south face [L ~> m]
  real, intent(in) :: dxCv_N         !< Cell x-length on north face [L ~> m]
  real, intent(in) :: dyCu_W         !< Cell y-length on west face [L ~> m]
  real, intent(in) :: dyCu_E         !< Cell y-length on east face [L ~> m]
  real, intent(in) :: rho            !< Ice density [R ~> kg m-3]
  real, intent(in) :: rhoi_rhow      !< Ice to ocean density ratio [nondim]
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

      dhdx_gp = ( (((-y_marginal_1) * h_nodal_cell(1,1)) + (( y_marginal_2) * h_nodal_cell(2,2))) + &
                  ((( y_marginal_1) * h_nodal_cell(2,1)) + ((-y_marginal_2) * h_nodal_cell(1,2))) ) / a
      dhdy_gp = ( (((-x_marginal_1) * h_nodal_cell(1,1)) + (( x_marginal_2) * h_nodal_cell(2,2))) + &
                  (((-x_marginal_2) * h_nodal_cell(2,1)) + (( x_marginal_1) * h_nodal_cell(1,2))) ) / d

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

      ! Local flotation test, even with GL_couple.
      is_grounded = (rhoi_rhow * h_gp - bed_gp > 0.0)

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

end subroutine calc_shelf_driving_stress_DG_subgrid

!> DG(1) driving-stress volume integral on the SEP2 partition of a grounding-line cell.
!! Each QP lies on one side of the cut and takes that side's grad(s) branch.
subroutine calc_shelf_driving_stress_DG_sep2(CS, h_nodal_cell, bed_corners, &
    dxCv_S, dxCv_N, dyCu_W, dyCu_E, rho, rhoi_rhow, grav, vol_dx, vol_dy, &
    sx_shelf, sy_shelf, calc_slope_diag)
  type(ice_shelf_dyn_CS), intent(in) :: CS    !< Ice shelf control structure
  real, dimension(2,2), intent(in) :: h_nodal_cell !< Q1 nodal thickness at the 4 corners [Z ~> m]
  real, dimension(2,2), intent(in) :: bed_corners  !< Bed depth at the 4 cell corners [Z ~> m]
  real, intent(in) :: dxCv_S         !< Cell x-length on south face [L ~> m]
  real, intent(in) :: dxCv_N         !< Cell x-length on north face [L ~> m]
  real, intent(in) :: dyCu_W         !< Cell y-length on west face [L ~> m]
  real, intent(in) :: dyCu_E         !< Cell y-length on east face [L ~> m]
  real, intent(in) :: rho            !< Ice density [R ~> kg m-3]
  real, intent(in) :: rhoi_rhow      !< rho/rhow [nondim]
  real, intent(in) :: grav           !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]
  real, dimension(2,2), intent(out) :: vol_dx !< Per-corner x volume integral [R L3 Z T-2 ~> kg m s-2]
  real, dimension(2,2), intent(out) :: vol_dy !< Per-corner y volume integral [R L3 Z T-2 ~> kg m s-2]
  real, intent(inout) :: sx_shelf    !< Cell-average x surface slope [Z L-1 ~> nondim]
  real, intent(inout) :: sy_shelf    !< Cell-average y surface slope [Z L-1 ~> nondim]
  logical, intent(in) :: calc_slope_diag !< True if slope diagnostics will be calculated

  real, dimension(4)     :: hc, bedc  ! Corner thickness and bed, flattened SW,SE,NW,NE [Z ~> m]
  real, dimension(4)     :: fls       ! Corner flotation deficit r*h - bed [Z ~> m]
  integer, dimension(4)  :: nqp       ! QPs per parent triangle
  real, dimension(4,7,4) :: beta      ! Corner-basis weights per (corner, QP, triangle) [nondim]
  real, dimension(7,4)   :: wref      ! Reference measure per (QP, triangle) [nondim]
  logical, dimension(7,4) :: qpg      ! Grounded state per (QP, triangle)
  real, dimension(4,7)   :: valx, valy ! Per-QP nodal contributions [R L3 Z T-2 ~> kg m s-2]
  real, dimension(4,4)   :: px, py    ! Per-(corner, triangle) partial sums [R L3 Z T-2 ~> kg m s-2]
  real, dimension(7) :: vsx, vsy      ! Per-QP Jacobian-weighted slopes [Z L ~> m2 m-1]
  real, dimension(7) :: vw            ! Per-QP Jacobian weights [L2 ~> m2]
  real, dimension(4) :: psx, psy      ! Per-triangle weighted-slope sums [Z L ~> m2 m-1]
  real, dimension(4) :: pw            ! Per-triangle weight sums [L2 ~> m2]
  real, dimension(4) :: dhxi, dheta   ! Per-triangle P1 gradient of h in reference space [Z ~> m]
  real, dimension(4) :: dbxi, dbeta   ! Per-triangle P1 gradient of bed in reference space [Z ~> m]
  real :: b1, b2, b3, b4    ! Corner-basis weights at the QP [nondim]
  real :: mS, mN, mW, mE    ! Marginal sums: interpolation weights of the 4 cell edges [nondim]
  real :: a, d              ! Interpolated cell-edge spacings at the QP [L ~> m]
  real :: weight            ! Quadrature weight wref * (a*d) [L2 ~> m2]
  real :: hloc              ! Ice thickness at the QP [Z ~> m]
  real :: dhdx_gp, dhdy_gp  ! Thickness gradients at the QP [Z L-1 ~> nondim]
  real :: dbdx_gp, dbdy_gp  ! Bed gradients at the QP [Z L-1 ~> nondim]
  real :: dsdx_gp, dsdy_gp  ! Surface gradients at the QP [Z L-1 ~> nondim]
  real :: fx_gp, fy_gp      ! Driving-stress integrand at the QP [R L Z T-2 ~> kg m-1 s-2]
  real :: w_total           ! Total Jacobian weight over the cell [L2 ~> m2]
  integer :: t, k, c

  hc(1) = h_nodal_cell(1,1) ; hc(2) = h_nodal_cell(2,1)
  hc(3) = h_nodal_cell(1,2) ; hc(4) = h_nodal_cell(2,2)
  bedc(1) = bed_corners(1,1) ; bedc(2) = bed_corners(2,1)
  bedc(3) = bed_corners(1,2) ; bedc(4) = bed_corners(2,2)

  fls(:) = (rhoi_rhow * hc(:)) - bedc(:)
  call sep2_cell_qps(fls, nqp, beta, wref, qpg)

  ! P1 fan gradients, the same interpolant sep2_cell_qps cut on.
  call sep2_fan_gradient(hc,   dhxi, dheta)
  call sep2_fan_gradient(bedc, dbxi, dbeta)

  do t=1,4
    do k=1,nqp(t)
      b1 = beta(1,k,t) ; b2 = beta(2,k,t) ; b3 = beta(3,k,t) ; b4 = beta(4,k,t)
      mS = b1 + b2 ; mN = b3 + b4 ; mW = b1 + b3 ; mE = b2 + b4
      a = (dxCv_S * mS) + (dxCv_N * mN)
      d = (dyCu_W * mW) + (dyCu_E * mE)
      weight = wref(k,t) * (a * d)

      hloc = ((b1 * hc(1)) + (b4 * hc(4))) + ((b2 * hc(2)) + (b3 * hc(3)))
      hloc = max(hloc, CS%min_h_shelf)
      dhdx_gp = dhxi(t) / a ; dhdy_gp = dheta(t) / d
      dbdx_gp = dbxi(t) / a ; dbdy_gp = dbeta(t) / d

      if (qpg(k,t)) then
        dsdx_gp = dhdx_gp - dbdx_gp
        dsdy_gp = dhdy_gp - dbdy_gp
      else
        dsdx_gp = (1.0 - rhoi_rhow) * dhdx_gp
        dsdy_gp = (1.0 - rhoi_rhow) * dhdy_gp
      endif

      fx_gp = -rho * grav * hloc * dsdx_gp
      fy_gp = -rho * grav * hloc * dsdy_gp

      do c=1,4
        valx(c,k) = (weight * beta(c,k,t)) * fx_gp
        valy(c,k) = (weight * beta(c,k,t)) * fy_gp
      enddo
      if (calc_slope_diag) then
        vw(k) = weight
        vsx(k) = dsdx_gp * weight
        vsy(k) = dsdy_gp * weight
      endif
    enddo

    ! Orbit-grouped QP sums: (2,3), (4,5) and (6,7) are reflection pairs.
    if (nqp(t) == 3) then
      do c=1,4
        px(c,t) = valx(c,1) + (valx(c,2) + valx(c,3))
        py(c,t) = valy(c,1) + (valy(c,2) + valy(c,3))
      enddo
      if (calc_slope_diag) then
        psx(t) = vsx(1) + (vsx(2) + vsx(3))
        psy(t) = vsy(1) + (vsy(2) + vsy(3))
        pw(t)  = vw(1) + (vw(2) + vw(3))
      endif
    else
      do c=1,4
        px(c,t) = (valx(c,1) + (valx(c,2) + valx(c,3))) + &
                  ((valx(c,4) + valx(c,5)) + (valx(c,6) + valx(c,7)))
        py(c,t) = (valy(c,1) + (valy(c,2) + valy(c,3))) + &
                  ((valy(c,4) + valy(c,5)) + (valy(c,6) + valy(c,7)))
      enddo
      if (calc_slope_diag) then
        psx(t) = (vsx(1) + (vsx(2) + vsx(3))) + ((vsx(4) + vsx(5)) + (vsx(6) + vsx(7)))
        psy(t) = (vsy(1) + (vsy(2) + vsy(3))) + ((vsy(4) + vsy(5)) + (vsy(6) + vsy(7)))
        pw(t)  = (vw(1) + (vw(2) + vw(3))) + ((vw(4) + vw(5)) + (vw(6) + vw(7)))
      endif
    endif
  enddo

  ! Role-grouped cross-triangle reduction (see CG_action_sep2_basal).
  vol_dx(1,1) = (px(1,1) + px(1,4)) + (px(1,2) + px(1,3)) ! SW: (S+W)+(E+N)
  vol_dx(2,1) = (px(2,2) + px(2,1)) + (px(2,3) + px(2,4)) ! SE: (E+S)+(N+W)
  vol_dx(1,2) = (px(3,4) + px(3,3)) + (px(3,1) + px(3,2)) ! NW: (W+N)+(S+E)
  vol_dx(2,2) = (px(4,3) + px(4,2)) + (px(4,4) + px(4,1)) ! NE: (N+E)+(W+S)
  vol_dy(1,1) = (py(1,1) + py(1,4)) + (py(1,2) + py(1,3))
  vol_dy(2,1) = (py(2,2) + py(2,1)) + (py(2,3) + py(2,4))
  vol_dy(1,2) = (py(3,4) + py(3,3)) + (py(3,1) + py(3,2))
  vol_dy(2,2) = (py(4,3) + py(4,2)) + (py(4,4) + py(4,1))

  if (calc_slope_diag) then
    ! Opposite-pair grouping is invariant under the triangle permutations of any
    ! rotation or reflection (S=1, E=2, N=3, W=4).
    w_total  = (pw(1) + pw(3)) + (pw(2) + pw(4))
    sx_shelf = ((psx(1) + psx(3)) + (psx(2) + psx(4))) / w_total
    sy_shelf = ((psy(1) + psy(3)) + (psy(2) + psy(4))) / w_total
  endif

end subroutine calc_shelf_driving_stress_DG_sep2

!> Set all DG(1) corners of cell (i,j) to h_shelf_value, for when h_shelf is overwritten
!! outside the advection.
subroutine reset_DG_to_cellmean_at_cell(CS, i, j, h_shelf_value)
  type(ice_shelf_dyn_CS), pointer    :: CS !< Ice shelf dynamics control structure.
  integer,                intent(in) :: i  !< i index of the cell to reset.
  integer,                intent(in) :: j  !< j index of the cell to reset.
  real,                   intent(in) :: h_shelf_value !< New cell-mean thickness to broadcast to all corners [Z ~> m]

  if (.not. associated(CS)) return
  if (.not. CS%use_DG_thickness) return
  CS%h_nodal(i,j,:,:) = h_shelf_value
end subroutine reset_DG_to_cellmean_at_cell

!> Pointer to CS%h_nodal, for INIT_ICE_THICKNESS_NODAL; null if CS is not associated.
function DG_nodal_thickness_ptr(CS) result(h_nodal)
  type(ice_shelf_dyn_CS),           pointer :: CS !< The ice shelf dynamics control structure
  real, dimension(:,:,:,:),         pointer :: h_nodal !< Q1 nodal thickness [Z ~> m]
  h_nodal => NULL()
  if (.not. associated(CS)) return
  h_nodal => CS%h_nodal
end function DG_nodal_thickness_ptr

!> Zero every DG(1) corner.
subroutine zero_DG_h_nodal(CS)
  type(ice_shelf_dyn_CS), pointer :: CS !< Ice shelf dynamics control structure.

  if (.not. associated(CS)) return
  if (.not. CS%use_DG_thickness) return
  CS%h_nodal(:,:,:,:) = 0.0
end subroutine zero_DG_h_nodal

!> Grounded area fraction of cell (i,j) from the partition the friction uses: the quadrant
!! fraction under CISM_FRICTION, otherwise CS%ground_frac (binary without a sub-element scheme).
pure function grounded_frac_cell(CS, i, j) result(fg)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< Ice shelf dynamics control structure.
  integer,                intent(in) :: i  !< i index of the cell.
  integer,                intent(in) :: j  !< j index of the cell.
  real :: fg !< The grounded area fraction of the cell [nondim]

  if (CS%gl_quad_friction) then
    fg = CS%f_ground_cell(i,j)
  else
    fg = CS%ground_frac(i,j)
  endif
  fg = min(max(fg, 0.0), 1.0)
end function grounded_frac_cell

!> Set ISS%water_flux from the MISMIP+ Ice1r melt profile (Leguy et al. 2021 eq. 18; Seroussi &
!! Morlighem 2018 eq. 4), m = min(max((d - 50 m)/15, 0), 30) m yr-1 with the flotation draft
!! d = (rho_i/rho_w)*h, scaled by CS%ice_only_melt_scale. Grounding-line cells follow
!! CS%ice_only_melt_glp. The thickness change itself is left to change_thickness_using_melt.
subroutine calc_prescribed_basal_melt(CS, ISS, G, US)
  type(ice_shelf_dyn_CS), pointer       :: CS  !< Ice shelf dynamics control structure.
  type(ice_shelf_state),  intent(inout) :: ISS !< Ice shelf state (hmask, h_shelf, water_flux).
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure.
  type(unit_scale_type),  intent(in)    :: US  !< A structure containing unit conversion factors

  real :: rhoi_rhow   ! The ratio of ice to ocean density [nondim]
  real :: draft       ! The flotation draft of the ice, (rho_i/rho_w)*h [Z ~> m]
  real :: melt_rate   ! The fully-floating basal melt rate, positive for melting [Z T-1 ~> m s-1]
  real :: fg          ! The grounded area fraction of the cell [nondim]
  real :: d_min       ! Draft at which melting begins, 50 m in Leguy eq. 18 [Z ~> m]
  real :: d_max       ! Draft at which melting saturates, 500 m in Leguy eq. 18 [Z ~> m]
  real :: m_max       ! The saturated melt rate, 30 m yr-1 in Leguy eq. 18 [Z T-1 ~> m s-1]
  real :: I_d_range   ! The reciprocal of (d_max - d_min) [Z-1 ~> m-1]
  logical :: floating ! True where the cell centre satisfies the flotation condition
  integer :: i, j

  if (.not. associated(CS)) return
  if (.not. CS%ice_only_basal_melt) return

  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  d_min = 50.0 * US%m_to_Z
  d_max = 500.0 * US%m_to_Z
  m_max = (30.0 / (365.0*86400.0)) * US%m_to_Z * US%T_to_s
  m_max = m_max * CS%ice_only_melt_scale
  I_d_range = 1.0 / (d_max - d_min)

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if ((ISS%hmask(i,j) /= 1.0) .and. (ISS%hmask(i,j) /= 2.0)) then
      ISS%water_flux(i,j) = 0.0 ; cycle
    endif

    draft = rhoi_rhow * ISS%h_shelf(i,j)
    melt_rate = m_max * min(max((draft - d_min) * I_d_range, 0.0), 1.0)

    !Apply according to chosen grounding line melt parameterization
    fg = grounded_frac_cell(CS, i, j)
    select case (CS%ice_only_melt_glp)
      case (MELT_GLP_FMP)
      case (MELT_GLP_FCMP)
        floating = (CS%bed_elev(i,j) - rhoi_rhow * max(ISS%h_shelf(i,j), CS%min_h_shelf)) >= 0.0
        if (.not. floating) melt_rate = 0.0
      case (MELT_GLP_PMP, MELT_GLP_SEM2)
        ! SEM2 has PMP's cell total; only its in-cell distribution differs.
        melt_rate = melt_rate * (1.0 - fg)
      case (MELT_GLP_NMP)
        if (fg > 0.0) melt_rate = 0.0
    end select

    ISS%water_flux(i,j) = melt_rate * CS%density_ice
  enddo ; enddo
end subroutine calc_prescribed_basal_melt

!> Add a cell-mean thickness source rate (positive for accumulation) at cell (i,j) to the
!! buffer applied by the next ice_shelf_advect_DG1_nodal, and to the basal buffer if basal.
!! Surface rate is recoverable as h_source_rate - h_source_rate_bmb
subroutine accumulate_DG_source_rate(CS, i, j, rate, basal)
  type(ice_shelf_dyn_CS), pointer    :: CS !< Ice shelf dynamics control structure.
  integer,                intent(in) :: i  !< i index of the cell.
  integer,                intent(in) :: j  !< j index of the cell.
  real,                   intent(in) :: rate !< Cell-mean thickness source rate
                                             !! to add [Z T-1 ~> m s-1].
  logical,                intent(in) :: basal !< If true, this rate is basal melt or
                                             !! freeze-on; if false it is surface mass balance.

  if (.not. associated(CS)) return
  if (.not. CS%use_DG_thickness) return
  CS%h_source_rate(i,j) = CS%h_source_rate(i,j) + rate
  if (basal) CS%h_source_rate_bmb(i,j) = CS%h_source_rate_bmb(i,j) + rate
end subroutine accumulate_DG_source_rate

!> Map a cell-mean source field to a Q1 nodal source with one SRC_OP_* operator:
!!  - SRC_OP_LOCAL: every corner of a cell takes that cell's own rate.
!!  - SRC_OP_AVERAGED: each B-grid corner takes the cell_mean_w-weighted average of its ice cells.
!!  - SRC_OP_SUBGRID: as AVERAGED but weighted by floating area (cell_mean_w*xi_basal), and each
!!    cell's corner gets that average times its own xi_basal.
subroutine project_source_to_nodes(CS, ISS, G, src, op, S_node, use_xi)
  type(ice_shelf_dyn_CS), intent(in)  :: CS  !< Ice shelf dynamics control structure.
  type(ice_shelf_state),  intent(in)  :: ISS !< Ice shelf state (hmask, h_shelf).
  type(ocean_grid_type),  intent(in)  :: G   !< The grid structure.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in) :: src !< Cell-mean source rate, halo
                                             !! updated by the caller [Z T-1 ~> m s-1].
  integer,                intent(in)  :: op  !< The cross-cell operator, SRC_OP_LOCAL,
                                             !! SRC_OP_AVERAGED or SRC_OP_SUBGRID.
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(out) :: S_node !< Q1 nodal source per
                                             !! cell at the 4 corners [Z T-1 ~> m s-1].
  logical,      optional, intent(in)  :: use_xi !< If true, distribute each cell's total within
                                             !! the cell in proportion to CS%xi_basal instead of
                                             !! uniformly (the SEM2 in-cell distribution).

  real, dimension(SZDIB_(G),SZDJB_(G)) :: S_corner ! Projected B-grid corner source [Z T-1 ~> m s-1]
  real, dimension(SZDIB_(G),SZDJB_(G)) :: w_corner ! Total cell_mean_w at a corner [L2 ~> m2]
  ! Slots by contributing cell (1=SW, 2=SE, 3=NW, 4=NE), summed in diagonal pairs for rotation invariance.
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: S_corner_b ! Per-cell corner source [Z L2 T-1 ~> m3 s-1]
  real, dimension(SZDIB_(G),SZDJB_(G),4) :: w_corner_b ! Per-cell corner weight [L2 ~> m2]
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: S_loc ! Nodal source before sharing [Z T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G)) :: m_eff ! Rate over the floating part of a cell [Z T-1 ~> m s-1]
  real :: w_contrib                                ! Corner contribution weight [L2 ~> m2]
  real :: w_sum                                    ! The DG area of a cell, sum_ab w [L2 ~> m2]
  real :: wxi_sum                                  ! The floating DG area, sum_ab w*xi [L2 ~> m2]
  logical :: xi_on                                 ! True when the SEM2 distribution is in use
  integer :: i, j, a, b

  xi_on = .false. ; if (present(use_xi)) xi_on = use_xi

  ! In-cell distribution: uniform, or with xi (proportional to nodal floating fraction)
  m_eff(:,:) = 0.0
  S_loc(:,:,:,:) = 0.0
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    if (ISS%hmask(i,j) /= 1.0) cycle
    if (xi_on) then
      w_sum = (CS%cell_mean_w(i,j,1,1) + CS%cell_mean_w(i,j,2,2)) + &
              (CS%cell_mean_w(i,j,2,1) + CS%cell_mean_w(i,j,1,2))
      wxi_sum = ((CS%cell_mean_w(i,j,1,1) * CS%xi_basal(i,j,1,1)) + &
                 (CS%cell_mean_w(i,j,2,2) * CS%xi_basal(i,j,2,2))) + &
                ((CS%cell_mean_w(i,j,2,1) * CS%xi_basal(i,j,2,1)) + &
                 (CS%cell_mean_w(i,j,1,2) * CS%xi_basal(i,j,1,2)))
      if (wxi_sum > 0.0) m_eff(i,j) = (w_sum * src(i,j)) / wxi_sum
      do b = 1, 2 ; do a = 1, 2
        S_loc(i,j,a,b) = m_eff(i,j) * CS%xi_basal(i,j,a,b)
      enddo ; enddo
    else
      m_eff(i,j) = src(i,j)
      S_loc(i,j,:,:) = src(i,j)
    endif
  enddo ; enddo

  S_node(:,:,:,:) = 0.0

  ! -- LOCAL --
  if (op == SRC_OP_LOCAL) then
    do j = G%jsc, G%jec ; do i = G%isc, G%iec
      if (ISS%hmask(i,j) /= 1.0) cycle
      S_node(i,j,:,:) = S_loc(i,j,:,:)
    enddo ; enddo
    return
  endif

  ! -- SUBGRID --
  S_corner(:,:) = 0.0
  w_corner(:,:) = 0.0
  S_corner_b(:,:,:) = 0.0
  w_corner_b(:,:,:) = 0.0

  if (op == SRC_OP_SUBGRID) then
    ! Average over the floating part of each corner's support; grounded corners (xi = 0) neither
    ! give nor receive. A uniform rate reduces to SRC_OP_LOCAL.
    do j = G%jsd, G%jed ; do i = G%isd, G%ied
      if (ISS%hmask(i,j) /= 1.0) cycle
      !SW corner
      w_contrib = CS%cell_mean_w(i,j,1,1) * CS%xi_basal(i,j,1,1)
      S_corner_b(I-1, J-1, 4) = w_contrib * m_eff(i,j)
      w_corner_b(I-1, J-1, 4) = w_contrib
      !SE corner
      w_contrib = CS%cell_mean_w(i,j,2,1) * CS%xi_basal(i,j,2,1)
      S_corner_b(I,   J-1, 3) = w_contrib * m_eff(i,j)
      w_corner_b(I,   J-1, 3) = w_contrib
      !NW corner
      w_contrib = CS%cell_mean_w(i,j,1,2) * CS%xi_basal(i,j,1,2)
      S_corner_b(I-1, J  , 2) = w_contrib * m_eff(i,j)
      w_corner_b(I-1, J  , 2) = w_contrib
      !NE corner
      w_contrib = CS%cell_mean_w(i,j,2,2) * CS%xi_basal(i,j,2,2)
      S_corner_b(I,   J  , 1) = w_contrib * m_eff(i,j)
      w_corner_b(I,   J  , 1) = w_contrib
    enddo ; enddo

    do j = G%JsdB, G%JedB ; do i = G%IsdB, G%IedB
      S_corner(I,J) = (S_corner_b(I,J,1) + S_corner_b(I,J,4)) + &
                      (S_corner_b(I,J,2) + S_corner_b(I,J,3))
      w_corner(I,J) = (w_corner_b(I,J,1) + w_corner_b(I,J,4)) + &
                      (w_corner_b(I,J,2) + w_corner_b(I,J,3))
      if (w_corner(I,J) > 0.0) S_corner(I,J) = S_corner(I,J) / w_corner(I,J)
    enddo ; enddo

    do j = G%jsc, G%jec ; do i = G%isc, G%iec
      if (ISS%hmask(i,j) /= 1.0) cycle
      S_node(i,j,1,1) = CS%xi_basal(i,j,1,1) * S_corner(I-1, J-1)
      S_node(i,j,2,1) = CS%xi_basal(i,j,2,1) * S_corner(I,   J-1)
      S_node(i,j,1,2) = CS%xi_basal(i,j,1,2) * S_corner(I-1, J  )
      S_node(i,j,2,2) = CS%xi_basal(i,j,2,2) * S_corner(I,   J  )
    enddo ; enddo
    return
  endif

  ! -- AVERAGED --
  ! Weighted contributions of each hmask=1 cell to its 4 B-grid corners.
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    if (ISS%hmask(i,j) /= 1.0) cycle
    !SW corner
    w_contrib = CS%cell_mean_w(i,j,1,1)
    S_corner_b(I-1, J-1, 4) = w_contrib * S_loc(i,j,1,1)
    w_corner_b(I-1, J-1, 4) = w_contrib
    !SE corner
    w_contrib = CS%cell_mean_w(i,j,2,1)
    S_corner_b(I,   J-1, 3) = w_contrib * S_loc(i,j,2,1)
    w_corner_b(I,   J-1, 3) = w_contrib
    !NW corner
    w_contrib = CS%cell_mean_w(i,j,1,2)
    S_corner_b(I-1, J  , 2) = w_contrib * S_loc(i,j,1,2)
    w_corner_b(I-1, J  , 2) = w_contrib
    !NE corner
    w_contrib = CS%cell_mean_w(i,j,2,2)
    S_corner_b(I,   J  , 1) = w_contrib * S_loc(i,j,2,2)
    w_corner_b(I,   J  , 1) = w_contrib
  enddo ; enddo

  ! Weighted average at each B-grid node; zero where no ice cell contributes.
  do j = G%JsdB, G%JedB ; do i = G%IsdB, G%IedB
    S_corner(I,J) = (S_corner_b(I,J,1) + S_corner_b(I,J,4)) + &
                    (S_corner_b(I,J,2) + S_corner_b(I,J,3))
    w_corner(I,J) = (w_corner_b(I,J,1) + w_corner_b(I,J,4)) + &
                    (w_corner_b(I,J,2) + w_corner_b(I,J,3))
    if (w_corner(I,J) > 0.0) S_corner(I,J) = S_corner(I,J) / w_corner(I,J)
  enddo ; enddo

  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle
    S_node(i,j,1,1) = S_corner(I-1, J-1)
    S_node(i,j,2,1) = S_corner(I,   J-1)
    S_node(i,j,1,2) = S_corner(I-1, J  )
    S_node(i,j,2,2) = S_corner(I,   J  )
  enddo ; enddo
end subroutine project_source_to_nodes

!> Build the nodal source S_node for ice_shelf_advect_DG1_nodal, projecting the basal and surface
!! parts with their own operators, or the combined buffer in one pass when their operators agree.
subroutine project_h_source_rate_to_nodes(CS, ISS, G, S_node)
  type(ice_shelf_dyn_CS), intent(in)    :: CS  !< Ice shelf dynamics control structure.
  type(ice_shelf_state),  intent(in)    :: ISS !< Ice shelf state (hmask, h_shelf).
  type(ocean_grid_type),  intent(in)    :: G   !< The grid structure.
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(out) :: S_node !< Q1 nodal source per
                                             !! cell at the 4 corners [Z T-1 ~> m s-1].

  real, dimension(SZDI_(G),SZDJ_(G)) :: src_smb ! Surface part of the source rate [Z T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: S_smb ! Surface part of the nodal source [Z T-1 ~> m s-1]
  integer :: surface_op ! The cross-cell operator applied to the surface source
  integer :: i, j, a, b

  surface_op = SRC_OP_AVERAGED
  if (CS%dg_surface_source_local) surface_op = SRC_OP_LOCAL

  call pass_var(CS%h_source_rate, G%domain)

  ! One pass when both parts share an operator and SEM2 is off.
  if ((CS%dg_basal_source_op == surface_op) .and. .not.CS%dg_basal_source_sem2) then
    call project_source_to_nodes(CS, ISS, G, CS%h_source_rate, surface_op, S_node)
    return
  endif

  ! Otherwise project the parts separately and sum. src_smb is formed on the data domain.
  call pass_var(CS%h_source_rate_bmb, G%domain)

  ! Calculate surface part
  src_smb(:,:) = 0.0
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    src_smb(i,j) = CS%h_source_rate(i,j) - CS%h_source_rate_bmb(i,j)
  enddo ; enddo

  call project_source_to_nodes(CS, ISS, G, CS%h_source_rate_bmb, CS%dg_basal_source_op, S_node, &
                               use_xi=CS%dg_basal_source_sem2)
  call project_source_to_nodes(CS, ISS, G, src_smb, surface_op, S_smb)

  do b = 1, 2 ; do a = 1, 2
    do j = G%jsc, G%jec ; do i = G%isc, G%iec
      S_node(i,j,a,b) = S_node(i,j,a,b) + S_smb(i,j,a,b)
    enddo ; enddo
  enddo ; enddo
end subroutine project_h_source_rate_to_nodes

!> Debug check that the global sum_ab cell_mean_w*S_node equals the global sum_ab cell_mean_w*S_cell,
!! i.e. that the source projection conserves mass.
subroutine check_nodal_source_conservation(CS, ISS, G, S_cell, S_node, label)
  type(ice_shelf_dyn_CS), intent(in) :: CS   !< Ice shelf dynamics control structure.
  type(ice_shelf_state),  intent(in) :: ISS  !< Ice shelf state (hmask).
  type(ocean_grid_type),  intent(in) :: G    !< The grid structure.
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in) :: S_cell !< The intended cell-mean source
                                             !! rate [Z T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(in) :: S_node !< The projected nodal source
                                             !! rate at the 4 corners [Z T-1 ~> m s-1].
  character(len=*),       intent(in) :: label !< Text identifying the caller in the message.

  real, dimension(SZDI_(G),SZDJ_(G)) :: tmp_node ! Per-cell nodal source integral [Z L2 T-1 ~> m3 s-1]
  real, dimension(SZDI_(G),SZDJ_(G)) :: tmp_cell ! Per-cell intended source integral [Z L2 T-1 ~> m3 s-1]
  real :: total_node ! Global integral of the projected nodal source [m3 s-1]
  real :: total_cell ! Global integral of the intended cell-mean source [m3 s-1]
  real :: denom      ! Larger of the two integrals in magnitude, for a relative error [m3 s-1]
  real :: w_sum      ! The DG area of a cell, sum_ab cell_mean_w [L2 ~> m2]
  real :: unscale    ! Conversion factor from [Z L2 T-1] to [m3 s-1]
  character(len=256) :: mesg
  integer :: i, j, is, ie, js, je, isr, ier, jsr, jer

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isr = is - (G%isd-1) ; ier = ie - (G%isd-1) ; jsr = js - (G%jsd-1) ; jer = je - (G%jsd-1)
  unscale = (G%US%Z_to_m * G%US%L_to_m**2) * G%US%s_to_T

  tmp_node(:,:) = 0.0 ; tmp_cell(:,:) = 0.0
  do j=js,je ; do i=is,ie
    if (ISS%hmask(i,j) /= 1.0) cycle
    tmp_node(i,j) = ((CS%cell_mean_w(i,j,1,1) * S_node(i,j,1,1)) + &
                     (CS%cell_mean_w(i,j,2,2) * S_node(i,j,2,2))) + &
                    ((CS%cell_mean_w(i,j,2,1) * S_node(i,j,2,1)) + &
                     (CS%cell_mean_w(i,j,1,2) * S_node(i,j,1,2)))
    w_sum = (CS%cell_mean_w(i,j,1,1) + CS%cell_mean_w(i,j,2,2)) + &
            (CS%cell_mean_w(i,j,2,1) + CS%cell_mean_w(i,j,1,2))
    tmp_cell(i,j) = w_sum * S_cell(i,j)
  enddo ; enddo

  total_node = reproducing_sum(tmp_node, isr, ier, jsr, jer, unscale=unscale)
  total_cell = reproducing_sum(tmp_cell, isr, ier, jsr, jer, unscale=unscale)

  denom = max(abs(total_cell), abs(total_node))
  if (is_root_pe()) then
    if (denom > 0.0) then
      write(mesg,'("DG source conservation (",A,"): nodal=",ES22.15," cell=",ES22.15, &
                  &" rel_err=",ES12.5)') trim(label), total_node, total_cell, &
                  (total_node - total_cell) / denom
    else
      write(mesg,'("DG source conservation (",A,"): both integrals are zero.")') trim(label)
    endif
    call MOM_mesg(trim(mesg))
  endif
end subroutine check_nodal_source_conservation

!> Debug check of sum_ab w*xi = (1 - ground_frac)*sum_ab w, which holds by partition of unity
!! (exactly on a uniform grid). A mis-indexed xi still conserves mass, so only this catches it.
subroutine check_xi_basal_consistency(CS, ISS, G)
  type(ice_shelf_dyn_CS), intent(in) :: CS   !< Ice shelf dynamics control structure.
  type(ice_shelf_state),  intent(in) :: ISS  !< Ice shelf state (hmask).
  type(ocean_grid_type),  intent(in) :: G    !< The grid structure.

  real :: w_sum      ! The DG area of a cell, sum_ab cell_mean_w [L2 ~> m2]
  real :: wxi_sum    ! The floating DG area, sum_ab cell_mean_w*xi [L2 ~> m2]
  real :: dev        ! Absolute deviation from the identity, as an area fraction [nondim]
  real :: dev_max    ! The largest deviation over the computational domain [nondim]
  character(len=256) :: mesg
  integer :: i, j

  dev_max = 0.0
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (ISS%hmask(i,j) /= 1.0) cycle
    w_sum = (CS%cell_mean_w(i,j,1,1) + CS%cell_mean_w(i,j,2,2)) + &
            (CS%cell_mean_w(i,j,2,1) + CS%cell_mean_w(i,j,1,2))
    if (w_sum <= 0.0) cycle
    wxi_sum = ((CS%cell_mean_w(i,j,1,1) * CS%xi_basal(i,j,1,1)) + &
               (CS%cell_mean_w(i,j,2,2) * CS%xi_basal(i,j,2,2))) + &
              ((CS%cell_mean_w(i,j,2,1) * CS%xi_basal(i,j,2,1)) + &
               (CS%cell_mean_w(i,j,1,2) * CS%xi_basal(i,j,1,2)))
    dev = abs((wxi_sum / w_sum) - (1.0 - CS%ground_frac(i,j)))
    dev_max = max(dev_max, dev)
  enddo ; enddo
  call max_across_PEs(dev_max)

  if (is_root_pe()) then
    write(mesg,'("Nodal floating fraction consistency: max |sum(w*xi)/sum(w) - (1-ground_frac)| =",&
                &ES12.5)') dev_max
    call MOM_mesg(trim(mesg))
  endif
end subroutine check_xi_basal_consistency


! ===========================================================================
! Nodal DG(1) helpers. CS%h_nodal is the DG thickness state; the per-cell
! metrics Minv_nodal and cell_mean_w are built in init_nodal_DG_metric.
! ===========================================================================

!> Read the DG(1) source, artificial-viscosity and mode-damping parameters.
subroutine read_DG_params(param_file, mdl, CS, US)
  type(param_file_type),   intent(in)    :: param_file
  character(len=*),        intent(in)    :: mdl
  type(ice_shelf_dyn_CS),  intent(inout) :: CS
  type(unit_scale_type),   intent(in)    :: US

  character(len=16) :: src_scheme_str ! DG(1) basal source cross-cell operator name


  call get_param(param_file, mdl, "DG_BASAL_SOURCE_SCHEME", src_scheme_str, &
                 "How the basal (melt) part of the DG(1) thickness source is shared between "//&
                 "cells at a corner. 'AVERAGED': each corner takes the weighted average of its "//&
                 "cells, which spreads melt across the grounding line. 'LOCAL': each corner takes "//&
                 "its own cell's rate. 'SUBGRID': the rate is averaged over the floating part of "//&
                 "each corner's support and returned by floating fraction, so grounded corners "//&
                 "neither give nor receive; requires ICE_ONLY_BASAL_MELT_GLP = 'SEM2'. All are "//&
                 "mass-conservative. Requires USE_DG_THICKNESS.", &
                 default="AVERAGED", do_not_log=.not.CS%use_DG_thickness)
  select case (trim(src_scheme_str))
    case ("AVERAGED") ; CS%dg_basal_source_op = SRC_OP_AVERAGED
    case ("LOCAL")    ; CS%dg_basal_source_op = SRC_OP_LOCAL
    case ("SUBGRID")  ; CS%dg_basal_source_op = SRC_OP_SUBGRID
    case default      ; call MOM_error(FATAL, "MOM_ice_shelf_dynamics: "//&
                          "DG_BASAL_SOURCE_SCHEME must be 'AVERAGED', 'LOCAL' or 'SUBGRID', "//&
                          "but got '"//trim(src_scheme_str)//"'.")
  end select
  if ((CS%dg_basal_source_op == SRC_OP_SUBGRID) .and. .not.CS%dg_basal_source_sem2) &
    call MOM_error(FATAL, "MOM_ice_shelf_dynamics: DG_BASAL_SOURCE_SCHEME = 'SUBGRID' requires "//&
                   "ICE_ONLY_BASAL_MELT_GLP = 'SEM2'; without it, it is identical to 'AVERAGED'.")

  call get_param(param_file, mdl, "DG_SURFACE_SOURCE_LOCAL", CS%dg_surface_source_local, &
                 "If true, each corner takes its own cell's surface mass balance. If false, "//&
                 "each corner takes the weighted average of its cells. Both are "//&
                 "mass-conservative. Requires USE_DG_THICKNESS.", &
                 default=.false., do_not_log=.not.CS%use_DG_thickness)


  call get_param(param_file, mdl, "DG1_ART_VISC_ADVECT_COEF", CS%dg_art_visc_advect_coef, &
                 "Coefficient on |u_face| in the DG(1) artificial-viscosity velocity "//&
                 "u_eff = ADVECT_COEF*|u_face| + STRAIN_COEF*eps_e_face*dx_perp.", &
                 units="nondim", default=0.0, &
                 do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG1_ART_VISC_STRAIN_COEF", CS%dg_art_visc_strain_coef, &
                 "Coefficient on eps_e_face*dx_perp in the DG(1) artificial-viscosity velocity, "//&
                 "where eps_e_face is the SSA effective strain rate at the face. It damps jumps "//&
                 "where |u_face| is small. 8*STRAIN_COEF = 1 sets the damping time to the "//&
                 "strain time.", &
                 units="nondim", default=0.0, &
                 do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG1_ART_VISC_ADVECT_L_REF", CS%dg_art_visc_advect_L_ref, &
                 "If positive, the |u_face| term of the DG(1) artificial viscosity becomes "//&
                 "ADVECT_COEF*|u_face|*dx_perp/L_REF, making its damping rate "//&
                 "resolution-independent. Non-positive uses ADVECT_COEF*|u_face|.", &
                 units="m", default=-1.0, scale=US%m_to_L, &
                 do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG1_ART_VISC_TAU_FLOOR", CS%dg_art_visc_tau_floor, &
                 "If positive, add dx_perp/TAU_FLOOR to the DG(1) artificial-viscosity u_eff, "//&
                 "giving every active face a jump decay rate of at least 8*C_MAX/TAU_FLOOR, "//&
                 "including on stagnant ice. Non-positive disables it.", &
                 units="s", default=-1.0, scale=US%s_to_T, &
                 do_not_log=.not.CS%use_DG_thickness)

  call get_param(param_file, mdl, "DG1_TILT_RELAX", CS%dg_tilt_relax, &
                 "If true, relax the DG(1) tilts and twist toward the ones the neighbouring "//&
                 "cell means imply, wherever the artificial viscosity acts, at a fraction of "//&
                 "its jump decay rate. The viscosity's jump penalty drives the thickness to a "//&
                 "continuous field with the same cell means, and such fields still admit an "//&
                 "alternating tilt that no face jump reveals; this gives the penalty a target "//&
                 "that excludes it. The reference is taken in surface form from neighbours of "//&
                 "the same flotation state, and a cell the grounding line crosses is left alone. "//&
                 "Moves no mass.", &
                 default=.false., do_not_log=.not.CS%use_DG_thickness)
  call get_param(param_file, mdl, "DG1_TILT_RELAX_FRAC", CS%dg_tilt_relax_frac, &
                 "The tilt-relaxation rate as a fraction of the artificial viscosity's jump "//&
                 "decay rate on the cell's faces along the same axis. Ignored unless "//&
                 "DG1_TILT_RELAX_U_CUT is zero.", &
                 units="nondim", default=1.0, do_not_log=.not.CS%dg_tilt_relax)
  call get_param(param_file, mdl, "DG1_TILT_RELAX_U_CUT", CS%dg_tilt_relax_u_cut, &
                 "If positive, the tilt-relaxation rate is the speed law max(U_CUT - "//&
                 "DG1_TILT_RELAX_U_CREDIT*|u_n|, 0)/dx along each axis, where u_n is the "//&
                 "cell-centred velocity component on that axis, and DG1_TILT_RELAX_FRAC is "//&
                 "ignored. A speed is the grid-invariant form, because both the rate the "//&
                 "transport supplies and the harm the mode does scale as 1/dx. The rate then "//&
                 "carries no dependence on the artificial viscosity, whose gate reads the face "//&
                 "jump and so is blind to this mode by construction, and is smallest where the "//&
                 "mode is purest. Ice moving faster than U_CUT/U_CREDIT along an axis is not "//&
                 "relaxed on that axis, because transport already removes the mode there.", &
                 units="m s-1", default=1.5854896E-06, scale=US%m_s_to_L_T, &
                 do_not_log=.not.CS%dg_tilt_relax)
  call get_param(param_file, mdl, "DG1_TILT_RELAX_U_CREDIT", CS%dg_tilt_relax_u_credit, &
                 "The fraction of the transport's own removal speed that DG1_TILT_RELAX_U_CUT "//&
                 "credits. Upwind transport removes the grid-scale tilt mode at |u_n|/dx, but "//&
                 "the artificial viscosity slows that to f(alpha)*|u_n|/dx, so a full credit "//&
                 "withdraws more than the flow returns. Under-crediting leaves the relaxation "//&
                 "on until |u_n| reaches U_CUT/U_CREDIT.", &
                 units="nondim", default=0.5, do_not_log=.not.CS%dg_tilt_relax)
  call get_param(param_file, mdl, "DG1_TILT_RELAX_TWIST", CS%dg_tilt_relax_twist, &
                 "If true, the tilt relaxation also relaxes the twist, toward the cross "//&
                 "difference of the four diagonal neighbours' means.", &
                 default=.true., do_not_log=.not.CS%dg_tilt_relax)
  call get_param(param_file, mdl, "DG1_TILT_RELAX_FREE_EDGE", CS%dg_tilt_relax_free_edge, &
                 "If true, the tilt relaxation also acts on an axis where the cell has no "//&
                 "artificial-viscosity face on one side, which happens at a true domain edge or "//&
                 "a calving front but not against a thickness boundary. Such a node has no "//&
                 "viscosity pin, so if the cell also has no usable reference on that axis nothing "//&
                 "holds it at all and it drifts to an unrealistically small thickness. The cell "//&
                 "is relaxed whatever its flotation state, using whichever reference the other "//&
                 "side gives: one of its own state if there is one, otherwise any ice cell. "//&
                 "Reading across a grounding line is allowed here, and only here, because the "//&
                 "alternative is no constraint at all.", &
                 default=.false., do_not_log=.not.CS%dg_tilt_relax)
  call get_param(param_file, mdl, "DG1_TILT_RELAX_WIDE_NEIGHBOURS", CS%dg_tilt_relax_wide_nb, &
                 "If true, a neighbour is usable for the tilt reference whenever it holds ice, "//&
                 "whatever its flotation state, and a cell the grounding line crosses is relaxed "//&
                 "like any other. Without this the relaxation does nothing along a continuous "//&
                 "curve the length of every grounding line, because a crossed cell is neither "//&
                 "relaxed nor usable as a neighbour. "//&
                 "a thickness reference cannot be read across a grounding line.", &
                 default=.false., do_not_log=.not.CS%dg_tilt_relax)

  call get_param(param_file, mdl, "DG1_TILT_DAMP", CS%dg_tilt_damp, &
                 "If true, damp the grid-scale DG(1) in-cell tilt. A tilt alternating between "//&
                 "cells gives no face jump, so the upwind flux and the artificial viscosity "//&
                 "cannot see it, yet it adds a spurious surface slope. The correction is zero "//&
                 "on a uniform tilt and moves no mass.", &
                 default=.true., do_not_log=(.not.CS%use_DG_thickness))
  if (.not.CS%use_DG_thickness) CS%dg_tilt_damp = .false.
  if (CS%use_DG_thickness .and. .not.CS%GL_regularize) call MOM_error(FATAL, &
    "USE_DG_THICKNESS requires GROUNDING_LINE_INTERPOLATE=True: the sub-element "//&
    "grounded fraction CS%ground_frac is what makes the grounding line sub-grid, "//&
    "and the mode damping grades itself on it.")

  call get_param(param_file, mdl, "DG1_TWIST_DAMP", CS%dg_twist_damp, &
                 "If true, also damp the grid-scale DG(1) in-cell twist, which is invisible "//&
                 "to the face jumps for the same reason as the tilt. It gives no net driving "//&
                 "stress. Uses DG1_TILT_DAMP_U_CUT and DG1_TILT_DAMP_R_HI.", &
                 default=.true., do_not_log=(.not.CS%use_DG_thickness))
  if (.not.CS%use_DG_thickness) CS%dg_twist_damp = .false.

  call get_param(param_file, mdl, "DG1_TILT_DAMP_ADVECTIVE", CS%dg_damp_advective, &
                 "If true, the mode damper supplies only the rate the transport does not "//&
                 "already deliver, per direction. Upwind transport removes the grid-scale "//&
                 "mode at about |u_n|/dx, so the damper mainly acts on slow ice.", &
                 default=.true., do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_ADVECTIVE_C", CS%dg_damp_advective_c, &
                 "Coefficient c in the mode-damper rates c*gate*U_CUT/dx and c*|u_n|/dx (the "//&
                 "latter credited by DG1_TILT_DAMP_ADVECTIVE). 1 is exact at the grid scale for "//&
                 "1D upwind advection; values below 1 credit the flow with less.", &
                 units="nondim", default=0.5, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_U_CUT", CS%dg_damp_u_cut, &
                 "Sets the mode damper's strength as a speed: the damping rate is the rate at "//&
                 "which upwind transport at speed gate*U_CUT would remove the mode, c*gate*U_CUT/dx. "//&
                 "With DG1_TILT_DAMP_ADVECTIVE the rate the flow already gives, c*|u_n|/dx, is "//&
                 "subtracted, so ice faster than gate*U_CUT is not damped. dx is the cell size "//&
                 "along the mode (sqrt(dx*dy) for the twist). Must be positive.", &
                 units="m s-1", default=6.341958E-06, scale=US%m_s_to_L_T, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_EXCESS_ONLY", CS%dg_damp_excess_only, &
                 "If true, the mode damper removes only the part of its detector that no "//&
                 "reference explains, instead of the whole detector once a threshold is "//&
                 "crossed. The two agree on a pure alternating mode.", &
                 default=.true., do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_KINK_REF", CS%dg_damp_kink_ref, &
                 "If true, the mode damper's reference reconstructs the slope break at the "//&
                 "grounding line from the grounded fraction, with rise f*s_g + (1-f)*s_f across "//&
                 "a cell, instead of exempting those cells. Where the branch slopes cannot be "//&
                 "built, it falls back to DG1_TILT_DAMP_GL_REACH.", &
                 default=.true., do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_KINK_TOL", CS%dg_damp_kink_tol, &
                 "Grounded-fraction misfit at which confidence in the twist's straight-line "//&
                 "reconstruction of the grounding line falls to zero, reverting to the "//&
                 "mean-supported reference.", &
                 units="nondim", default=0.2, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_KINK_FIT_TOL", CS%dg_damp_kink_fit_tol, &
                 "Relative slope misfit at which confidence in the tilt's two-branch "//&
                 "reconstruction of the grounding line falls to zero, reverting the reference "//&
                 "and closing the gate. Non-positive skips the test.", &
                 units="nondim", default=-1.0, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_GL_REACH", CS%dg_damp_gl_reach, &
                 "Reach of the mode damper's grounding-line protection: 0 the cell the "//&
                 "flotation contour crosses, 1 the cells in the detector stencil, 2 the full "//&
                 "reference stencil. 1 and 2 test whether the window has both grounded and "//&
                 "floating ice.", &
                 default=1, do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))
  if ((CS%dg_damp_gl_reach < 0) .or. (CS%dg_damp_gl_reach > 2)) call MOM_error(FATAL, &
    "DG1_TILT_DAMP_GL_REACH must be 0, 1 or 2.")


  if ((CS%dg_tilt_damp .or. CS%dg_twist_damp) .and. (CS%dg_damp_u_cut <= 0.0)) &
    call MOM_error(FATAL, "MOM_ice_shelf_dynamics: DG1_TILT_DAMP_U_CUT must be positive; "//&
                   "it sets the damping rate, which would otherwise be zero.")

  call get_param(param_file, mdl, "DG1_TILT_DAMP_DETECTOR", CS%dg_damp_detector, &
                 "Which detector finds the grid-scale zigzag.  0 is the tilt Laplacian held "//&
                 "against a bed floor, a mean-supported floor and a grounding-line kink "//&
                 "reference.  1 is stencil agreement: read the second difference of the slope "//&
                 "on each group of three adjacent cells, on the thickness and on the surface "//&
                 "form, and keep only the part every reading agrees on.  1 has no parameters "//&
                 "and gives exactly zero for any structure four cells wide or wider.  2 lets "//&
                 "each cell choose: 0 where the cells the stencil reads are all grounded or "//&
                 "all afloat, and 1 where they are not.  0 finds more of the mode, but it "//&
                 "needs a reference that averages across a flotation break, which is what "//&
                 "the kink reference exists to repair; 1 needs no reference at all.  2 "//&
                 "therefore uses 0 where it is strong and 1 where it is safe, and the kink "//&
                 "reference can be left off.", &
                 default=DAMP_DET_TODAY, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))
  if ((CS%dg_damp_detector < DAMP_DET_TODAY) .or. &
      (CS%dg_damp_detector > DAMP_DET_HYBRID)) call MOM_error(FATAL, &
    "DG1_TILT_DAMP_DETECTOR must be 0, 1 or 2.")

  call get_param(param_file, mdl, "DG1_TILT_DAMP_EDGE", CS%dg_damp_edge_rule, &
                 "What the stencil-agreement detector does at a cell that has lost a "//&
                 "neighbour, next to an ice front, a wall or a nunatak. The cancellation on "//&
                 "smooth ice works through the opposite signs of the centred reading and the "//&
                 "one-sided ones, so it needs a neighbour on both sides; a cell without one "//&
                 "also has no spread, so its gate cannot close. 0 leaves those cells to the "//&
                 "readings they have, which damps smooth ice at the full rate. 1 takes the "//&
                 "estimates from the cell means instead, over every run of 2 to 5 usable "//&
                 "cells that holds the cell, and stands down where that run crosses a "//&
                 "flotation break or holds fewer than 3 cells. The reference is exact for "//&
                 "cells of equal width and first order otherwise.", &
                 default=DAMP_EDGE_OFF, do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)) &
                                  .or. ((CS%dg_damp_detector /= DAMP_DET_AGREE) .and. &
                                       (CS%dg_damp_detector /= DAMP_DET_HYBRID)))
  if ((CS%dg_damp_edge_rule < 0) .or. (CS%dg_damp_edge_rule > 2)) call MOM_error(FATAL, &
    "read_DG_params: DG1_TILT_DAMP_EDGE must be 0, 1 or 2.")
  call get_param(param_file, mdl, "DG1_TILT_DAMP_FIELDS", CS%dg_damp_fields, &
                 "Which fields the stencil-agreement detector reads. 0 reads the thickness "//&
                 "and the surface form, and acts only where both agree; a thickness that "//&
                 "follows a rough bed is then protected, because the surface does not see "//&
                 "the bed. 1 reads the surface form alone. 2 reads the thickness alone. "//&
                 "Where the ice floats the two are the same field, so this changes grounded "//&
                 "ice only. Setting 0 is the careful one, and it is what stops the damper "//&
                 "eating a thickness that follows a rough bed. 3 keeps that everywhere "//&
                 "except at a cell the flotation contour crosses: there the surface form is "//&
                 "a FRACTION of the bed removed, which is no surface at all, so it is "//&
                 "dropped and the thickness decides alone.", &
                 default=DAMP_FLD_BOTH, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)) &
                            .or. ((CS%dg_damp_detector /= DAMP_DET_AGREE) .and. &
                                  (CS%dg_damp_detector /= DAMP_DET_HYBRID)))
  if ((CS%dg_damp_fields < DAMP_FLD_BOTH) .or. (CS%dg_damp_fields > DAMP_FLD_GL_TH)) &
    call MOM_error(FATAL, "read_DG_params: DG1_TILT_DAMP_FIELDS must be 0, 1, 2 or 3.")
  call get_param(param_file, mdl, "DG1_TILT_DAMP_CENTRED_AMP", CS%dg_damp_centred_amp, &
                 "If true, the minmod of the stencil-agreement readings only decides whether "//&
                 "to damp, and the centred reading sets how much is removed. The minmod keeps "//&
                 "the veto, so a structure 4 cells wide and wider still gives exactly zero, "//&
                 "and the response to the grid-scale mode itself does not change. Only the "//&
                 "band between 2 and 4 cells gains. If false, the minmod sets the amount as "//&
                 "well, which makes the removal smaller than the mode whenever the stencils "//&
                 "read different mixtures of the mode and of real structure.", &
                 default=.false., do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)) &
                                  .or. ((CS%dg_damp_detector /= DAMP_DET_AGREE) .and. &
                                       (CS%dg_damp_detector /= DAMP_DET_HYBRID)))
  call get_param(param_file, mdl, "DG1_TILT_DAMP_REDUCE", CS%dg_damp_reduce, &
                 "How the stencil-agreement detector reduces its readings to one value. 0 "//&
                 "compares every reading with every other one, so a single reading of the "//&
                 "opposite sign stops the damper. The damper then drives the ice until one "//&
                 "reading sits on zero and stops itself. 1 adds the two one-sided readings "//&
                 "before it compares them with the centred one, which ends that, because the "//&
                 "sum of the two does not reach zero where one of them does. On smooth ice "//&
                 "the sum still takes the sign opposite to the centred reading, so a "//&
                 "structure 4 cells wide and wider still gives exactly zero. The pair is "//&
                 "used only where the cells the stencil reads share a grounding status: one "//&
                 "cell of grounded ice beside floating ice puts a real kink in the surface "//&
                 "that reads exactly like the mode, and there setting 0 is what leaves the "//&
                 "grounding line in place.", &
                 default=DAMP_RED_MINMOD, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)) &
                            .or. ((CS%dg_damp_detector /= DAMP_DET_AGREE) .and. &
                                       (CS%dg_damp_detector /= DAMP_DET_HYBRID)))
  if ((CS%dg_damp_reduce < DAMP_RED_MINMOD) .or. (CS%dg_damp_reduce > DAMP_RED_PAIRED)) &
    call MOM_error(FATAL, "read_DG_params: DG1_TILT_DAMP_REDUCE must be 0 or 1.")
  call get_param(param_file, mdl, "DG1_TILT_DAMP_SINGLE_RULE", CS%dg_damp_single_rule, &
                 "If true, the agreement detector leaves a tilt alone in a cell that has one "//&
                 "stencil, when that stencil is neither wholly floating nor wholly grounded.  "//&
                 "Such a cell is at an ice front or a wall and has no second stencil to "//&
                 "compare with, so a stencil that crosses the grounding line reads the change "//&
                 "of slope there.  The twist does not need this: it reads both axes.", &
                 default=.true., &
                 do_not_log=((CS%dg_damp_detector /= DAMP_DET_AGREE) .and. &
                             (CS%dg_damp_detector /= DAMP_DET_HYBRID)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_FILTER", CS%dg_damp_filter, &
                 "If true, apply a 1-2-1 high pass to the damper's removal before it is "//&
                 "applied, along the direction of a tilt and along both axes for the twist.  "//&
                 "The removal is an alternating pattern times a slowly changing envelope, and "//&
                 "the slowly changing part moves the smooth field; the high pass takes that "//&
                 "part out and leaves a pure zigzag exactly as it is.  It widens the damper's "//&
                 "reach by one cell.", &
                 default=.false., &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_GATE", CS%dg_damp_gate_form, &
                 "Which gate scales the mode damping.  0 compares the detected zigzag with the "//&
                 "ice thickness, through DG1_TILT_DAMP_R_HI; that gate closes on thick ice, and "//&
                 "the slope error it accepts doubles each time the grid is refined by two.  1 "//&
                 "compares it with the real surface slope, through DG1_TILT_DAMP_RHO_S; that "//&
                 "gate closes on steep grounded ice.  2 compares the part the detector readings "//&
                 "agree on with the part they do not, through DG1_TILT_DAMP_RHO_G; it carries "//&
                 "no thickness, no slope and no grid spacing, and it needs the agreement "//&
                 "detector to be meaningful.", &
                 default=DAMP_GATE_THICKNESS, &
                 do_not_log=(.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)))
  if ((CS%dg_damp_gate_form < DAMP_GATE_THICKNESS) .or. &
      (CS%dg_damp_gate_form > DAMP_GATE_AGREEMENT)) call MOM_error(FATAL, &
    "DG1_TILT_DAMP_GATE must be 0, 1 or 2.")
  ! Under the hybrid this names the gate of the tilt-Laplacian branch alone; the agreement
  ! branch always uses the agreement gate, which is the only one that reads its spread.  The
  ! agreement gate on the tilt Laplacian would stand fully open, because that detector reports
  ! no spread.
  if ((CS%dg_damp_detector == DAMP_DET_HYBRID) .and. &
      (CS%dg_damp_gate_form == DAMP_GATE_AGREEMENT)) call MOM_error(FATAL, &
    "read_DG_params: DG1_TILT_DAMP_GATE = 2 cannot gate the tilt-Laplacian branch of "//&
    "DG1_TILT_DAMP_DETECTOR = 2, which reports no spread.  Use 0 or 1; the agreement "//&
    "branch gates itself.")

  call get_param(param_file, mdl, "DG1_TILT_DAMP_RHO_S", CS%dg_damp_rho_s, &
                 "Zigzag surface slope, as a multiple of the real surface slope, at which "//&
                 "DG1_TILT_DAMP_GATE = 1 is fully open.  The detector reads twice the amplitude "//&
                 "of the zigzag, so 2 opens the gate where the zigzag is as steep as the ice.", &
                 units="nondim", default=2.0, &
                 do_not_log=(CS%dg_damp_gate_form /= DAMP_GATE_SLOPE))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_SLOPE_FLOOR", CS%dg_damp_slope_floor, &
                 "Smallest real surface slope that DG1_TILT_DAMP_GATE = 1 accepts as a "//&
                 "reference.  It is a slope rather than a rise, so the gate carries no grid "//&
                 "spacing on ice flat enough for the floor to decide.", &
                 units="nondim", default=1.0e-5, &
                 do_not_log=(CS%dg_damp_gate_form /= DAMP_GATE_SLOPE))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_RHO_G", CS%dg_damp_rho_g, &
                 "Disagreement between the detector readings, as a multiple of the part they "//&
                 "agree on, at which DG1_TILT_DAMP_GATE = 2 is half open.  Raising it protects "//&
                 "real features and costs almost no zigzag removal; lowering it does not "//&
                 "remove more zigzag, because what stays comes from the detector.", &
                 units="nondim", default=4.0, &
                 do_not_log=(CS%dg_damp_gate_form /= DAMP_GATE_AGREEMENT))

  call get_param(param_file, mdl, "DG1_TILT_DAMP_R_HI", CS%dg_tilt_damp_r_hi, &
                 "Gate value at which the mode damping is at full strength. The gate is "//&
                 "(ds/dh)*(|A| - max(|A_bed|,|A_ref|))/h, where A is the discrete Laplacian of "//&
                 "the cell tilt, A_bed and A_ref the same for the bed and for the tilt implied "//&
                 "by the cell means, and ds/dh is 1 grounded and 1 - rho_i/rho_w floating.  "//&
                 "Used only by DG1_TILT_DAMP_GATE = 0.", &
                 units="nondim", default=0.02, &
                 do_not_log=(CS%dg_damp_gate_form /= DAMP_GATE_THICKNESS))

  call get_param(param_file, mdl, "DG1_ART_VISC_R_HI", CS%dg_art_visc_r_hi, &
                 "Face surface jump, relative to the mean thickness, at which the DG(1) "//&
                 "artificial viscosity reaches its full coefficient.", &
                 units="nondim", default=0.005, &
                 do_not_log=.not.CS%use_DG_thickness)
  call get_param(param_file, mdl, "DG1_ART_VISC_KCELL", CS%dg_art_visc_kcell, &
                 "Bound on dt times a cell's summed face jump decay rates "//&
                 "(4*amp*c*u_eff/dx_perp); the face coefficients are scaled down to meet it. "//&
                 "SSP-RK2 requires < 2.", &
                 units="nondim", default=1.0, &
                 do_not_log=.not.CS%use_DG_thickness)
  if (CS%use_DG_thickness .and. (CS%dg_art_visc_r_hi <= 0.0)) call MOM_error(FATAL, &
      "MOM_ice_shelf_dynamics, initialize_ice_shelf_dyn: DG1_ART_VISC_R_HI must be positive.")

  ! Stagnant-jump diagnostic thresholds.
  CS%dg_slow_idle_u_tiny   = 1.0e-9  * US%m_s_to_L_T
  CS%dg_slow_idle_eps_tiny = 1.0e-12 * US%s_to_T
  CS%dg_slow_idle_s_tol    = 1.0     * US%m_to_Z

end subroutine read_DG_params

!> Build the per-cell tables Minv_nodal and cell_mean_w. On an orthogonal grid the Q1 mass
!! matrix is M_xi (x) M_eta, with M_xi = int N_a N_a' d(xi) dxi, M_eta = int N_b N_b' a(eta) deta,
!! d(xi) = dyW*(1-xi) + dyE*xi and a(eta) = dxS*(1-eta) + dxN*eta. Its inverse is the product of
!! the analytic 2x2 inverses.
subroutine init_nodal_DG_metric(CS, G)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ocean_grid_type),  intent(in)    :: G

  real :: dxS, dxN, dyW, dyE    ! Face lengths [L ~> m]
  real :: M11, M12, M22         ! 1D mass-matrix entries [L ~> m]
  real :: det                   ! 1D mass-matrix determinant [L2 ~> m2]
  real, dimension(2,2) :: Minv_xi  ! Inverse 1D mass matrix along xi [L-1 ~> m-1]
  real, dimension(2,2) :: Minv_eta ! Inverse 1D mass matrix along eta [L-1 ~> m-1]
  integer :: i, j, isd, ied, jsd, jed, a, b, p, q

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  do j = jsd, jed ; do i = isd, ied
    call cell_face_lengths(G, i, j, CS%reentrant_x, CS%reentrant_y, dxS, dxN, dyW, dyE)

    !M_xi: int_0^1 N_a*N_a'*d(xi) dxi where d(xi) = dyW*(1-xi) + dyE*xi.
    ! Closed-form: M11 = dyW/4 + dyE/12, M22 = dyW/12 + dyE/4, M12 = dyW/12 + dyE/12.
    M11 = (dyW/4.0) + (dyE/12.0)
    M22 = (dyW/12.0) + (dyE/4.0)
    M12 = (dyW + dyE)/12.0
    det = M11*M22 - M12*M12
    Minv_xi(:,:) = 0.0
    if (det > 0.0) then
      Minv_xi(1,1) =  M22 / det
      Minv_xi(2,2) =  M11 / det
      Minv_xi(1,2) = -M12 / det
      Minv_xi(2,1) = -M12 / det
    endif

    ! M_eta: int_0^1 N_b*N_b'*a(eta) deta where a(eta) = dxS*(1-eta) + dxN*eta.
    M11 = (dxS/4.0) + (dxN/12.0)
    M22 = (dxS/12.0) + (dxN/4.0)
    M12 = (dxS + dxN)/12.0
    det = M11*M22 - M12*M12
    Minv_eta(:,:) = 0.0
    if (det > 0.0) then
      Minv_eta(1,1) =  M22 / det
      Minv_eta(2,2) =  M11 / det
      Minv_eta(1,2) = -M12 / det
      Minv_eta(2,1) = -M12 / det
    endif

    do q = 1, 2 ; do p = 1, 2 ; do b = 1, 2 ; do a = 1, 2
      CS%Minv_nodal(i,j,a,b,p,q) = Minv_xi(a,p) * Minv_eta(b,q)
    enddo ; enddo ; enddo ; enddo

    call corner_cell_weights(dxS, dxN, dyW, dyE, CS%cell_mean_w(i,j,:,:))
  enddo ; enddo
end subroutine init_nodal_DG_metric

!> At initialization on a symmetric reentrant grid, copy the DG corners across the wrap so both
!! cells at a periodic boundary node start with the same value.
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

  ! Collective, so every PE calls it.
  call pass_corner_field(CS%h_nodal, G)
end subroutine enforce_wrap_corner_consistency

!> Halo update of a per-cell corner field, each corner as a cell-centered scalar.
subroutine pass_corner_field(h_nodal, G)
  type(ocean_grid_type),  intent(in) :: G
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(inout) :: h_nodal
  real, dimension(SZDI_(G),SZDJ_(G),4) :: tmp
  real, dimension(SZDI_(G),SZDJ_(G)) :: hbar ! Plain corner mean of the cell [A ~> a]
  real, dimension(SZDI_(G),SZDJ_(G)) :: t_xi, t_eta ! Corner tilts along each axis [A ~> a]
  real, dimension(SZDI_(G),SZDJ_(G)) :: t_w ! Corner twist of the cell [A ~> a]
  integer :: i, j

  if (G%Domain%Y_FLAGS /= FOLD_NORTH_EDGE) then
    ! No fold: the four corners are four scalars at the cell centre and pass_var moves them
    ! correctly.  This path is kept because it is exact, where the decomposition below is not.
    tmp(:,:,1) = h_nodal(:,:,1,1)
    tmp(:,:,2) = h_nodal(:,:,2,1)
    tmp(:,:,3) = h_nodal(:,:,1,2)
    tmp(:,:,4) = h_nodal(:,:,2,2)
    call pass_var(tmp, G%domain)
    h_nodal(:,:,1,1) = tmp(:,:,1)
    h_nodal(:,:,2,1) = tmp(:,:,2)
    h_nodal(:,:,1,2) = tmp(:,:,3)
    h_nodal(:,:,2,2) = tmp(:,:,4)
    return
  endif

  ! Across a tripolar fold the halo cell (i, j) is filled from (isg+ieg-i, 2*jeg-j+1), so BOTH
  ! local axes reverse and its frame is turned through half a revolution: its south-west corner is
  ! the source cell's north-east corner.  pass_var gets the cell-to-cell mapping right but cannot
  ! permute inside a cell, because nothing tells it the four fields belong together, so the corner
  ! labels come back wrong and both tilts change sign.  Send the parts instead, each in the form
  ! whose behaviour at the fold the domain layer already knows:
  !   the mean and the twist do not change under a half turn, so they travel as scalars;
  !   the two tilts are a gradient, so they change sign together and travel as a vector pair.
  ! The half turn is a rotation and not a reflection -- the two axis reversals compose -- which is
  ! what lets the twist travel as a scalar.  A reflection would change its sign.
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    hbar(i,j)  = 0.25 * ((h_nodal(i,j,1,1) + h_nodal(i,j,2,2)) + &
                         (h_nodal(i,j,2,1) + h_nodal(i,j,1,2)))
    t_xi(i,j)  = 0.5 * ((h_nodal(i,j,2,1) - h_nodal(i,j,1,1)) + &
                        (h_nodal(i,j,2,2) - h_nodal(i,j,1,2)))
    t_eta(i,j) = 0.5 * ((h_nodal(i,j,1,2) - h_nodal(i,j,1,1)) + &
                        (h_nodal(i,j,2,2) - h_nodal(i,j,2,1)))
    t_w(i,j)   = (h_nodal(i,j,2,2) - h_nodal(i,j,1,2)) - &
                 (h_nodal(i,j,2,1) - h_nodal(i,j,1,1))
  enddo ; enddo

  call pass_var(hbar, G%domain, complete=.false.)
  call pass_var(t_w, G%domain, complete=.true.)
  call pass_vector(t_xi, t_eta, G%domain, stagger=AGRID)

  ! Rebuild.  This inverts the four lines above exactly.
  do j = G%jsd, G%jed ; do i = G%isd, G%ied
    h_nodal(i,j,1,1) = ((hbar(i,j) - (0.5*t_xi(i,j))) - (0.5*t_eta(i,j))) + (0.25*t_w(i,j))
    h_nodal(i,j,2,1) = ((hbar(i,j) + (0.5*t_xi(i,j))) - (0.5*t_eta(i,j))) - (0.25*t_w(i,j))
    h_nodal(i,j,1,2) = ((hbar(i,j) - (0.5*t_xi(i,j))) + (0.5*t_eta(i,j))) - (0.25*t_w(i,j))
    h_nodal(i,j,2,2) = ((hbar(i,j) + (0.5*t_xi(i,j))) + (0.5*t_eta(i,j))) + (0.25*t_w(i,j))
  enddo ; enddo

end subroutine pass_corner_field


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


!> Positivity limiter: scale each cell's corner deviations from its mean so that the smallest
!! corner is at least CS%min_h_shelf. The cell mean is unchanged.
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

!> (p1 + p2) + (p3 + p4) with each product rounded before the sum. A fused multiply-add keeps
!! one product unrounded, and a quarter turn cycles the corners through the slots, so no
!! ordering is rotation-invariant. Parentheses do not stop gfortran fusing once it vectorizes;
!! the volatile store does. Without FMA this equals the plain expression.
function dg_sum4_rounded(p1, p2, p3, p4) result(s)
  real, intent(in) :: p1 !< First product, paired with p2 [arbitrary]
  real, intent(in) :: p2 !< Second product, paired with p1 [arbitrary]
  real, intent(in) :: p3 !< Third product, paired with p4 [arbitrary]
  real, intent(in) :: p4 !< Fourth product, paired with p3 [arbitrary]
  real :: s              !< (p1 + p2) + (p3 + p4), each term rounded first [arbitrary]
  real, volatile, dimension(4) :: t ! The products, rounded by storage [arbitrary]
  t(1) = p1 ; t(2) = p2 ; t(3) = p3 ; t(4) = p4
  s = (t(1) + t(2)) + (t(3) + t(4))
end function dg_sum4_rounded

!> Apply a cell's inverse mass matrix, out(a,b) = sum_pq Minv_cell(a,b,p,q)*rhs(p,q).
subroutine apply_nodal_DG_mass_inverse(Minv_cell, rhs, out)
  real, dimension(2,2,2,2), intent(in)  :: Minv_cell !< Inverse mass matrix of this cell [L-2 ~> m-2]
  real, dimension(2,2),     intent(in)  :: rhs       !< Per-cell RHS at the 4 corners [Z L2 T-1 ~> m3 s-1]
  real, dimension(2,2),     intent(out) :: out       !< M^-1 * rhs [Z T-1 ~> m s-1]
  integer :: a, b
  ! Opposite corners paired for rotation invariance.
  do b = 1, 2 ; do a = 1, 2
    out(a,b) = dg_sum4_rounded(Minv_cell(a,b,1,1) * rhs(1,1), Minv_cell(a,b,2,2) * rhs(2,2), &
                               Minv_cell(a,b,2,1) * rhs(2,1), Minv_cell(a,b,1,2) * rhs(1,2))
  enddo ; enddo
end subroutine apply_nodal_DG_mass_inverse

!> Surface jump s_B - s_A at a face point, each side taking its own flotation branch.
pure function dg1_wb_surface_jump(h_A, h_B, bed_qp, rhoi_rhow) result(ds)
  real, intent(in) :: h_A       !< Side-A face-QP thickness [Z ~> m]
  real, intent(in) :: h_B       !< Side-B face-QP thickness [Z ~> m]
  real, intent(in) :: bed_qp    !< Bed elevation at the face QP, single-valued [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio [nondim]
  real :: ds                    !< Surface-elevation jump s_B - s_A [Z ~> m]
  real :: s_A, s_B              ! Per-side surface elevation [Z ~> m]
  real :: one_m_r               ! 1 - rhoi_rhow [nondim]
  one_m_r = 1.0 - rhoi_rhow
  if (rhoi_rhow*h_A - bed_qp > 0.0) then ; s_A = (h_A - bed_qp)
  else ; s_A = (one_m_r*h_A) ; endif
  if (rhoi_rhow*h_B - bed_qp > 0.0) then ; s_B = (h_B - bed_qp)
  else ; s_B = (one_m_r*h_B) ; endif
  ds = s_B - s_A
end function dg1_wb_surface_jump

!> Harmonic mean across a face of dh/ds (1 grounded, 1/(1-r) floating), which converts a surface
!! jump to a thickness jump and makes DG1_WB_JUMP_RATE_AMP exactly 2.
pure function dg1_wb_slope_mean(h_A, h_B, bed_qp, rhoi_rhow) result(m)
  real, intent(in) :: h_A       !< Side-A face-QP thickness [Z ~> m]
  real, intent(in) :: h_B       !< Side-B face-QP thickness [Z ~> m]
  real, intent(in) :: bed_qp    !< Bed elevation at the face QP, single-valued [Z ~> m]
  real, intent(in) :: rhoi_rhow !< Ice/ocean density ratio [nondim]
  real :: m                     !< Harmonic mean inverse flotation slope dh/ds [nondim]
  real :: g_A, g_B              ! Per-side surface slope ds/dh [nondim]
  real :: one_m_r               ! 1 - rhoi_rhow [nondim]
  one_m_r = 1.0 - rhoi_rhow
  if (rhoi_rhow*h_A - bed_qp > 0.0) then ; g_A = 1.0 ; else ; g_A = one_m_r ; endif
  if (rhoi_rhow*h_B - bed_qp > 0.0) then ; g_B = 1.0 ; else ; g_B = one_m_r ; endif
  if (g_A == g_B) then
    m = 1.0/g_A
  else
    m = 2.0/(g_A + g_B)
  endif
end function dg1_wb_slope_mean

!> SSA effective strain rate sqrt(eps_xx^2 + eps_yy^2 + eps_xx*eps_yy + eps_xy^2).
!! For DG(1) artificial-viscosity strain-scaled floor.
pure function dg1_face_eps_eff(dudx, dudy, dvdx, dvdy) result(eps_e)
  real, intent(in) :: dudx     !< du/dx at face midpoint [T-1 ~> s-1]
  real, intent(in) :: dudy     !< du/dy at face midpoint [T-1 ~> s-1]
  real, intent(in) :: dvdx     !< dv/dx at face midpoint [T-1 ~> s-1]
  real, intent(in) :: dvdy     !< dv/dy at face midpoint [T-1 ~> s-1]
  real :: eps_e                !< Effective strain rate [T-1 ~> s-1]
  real :: eps_xy               ! Shear strain rate [T-1 ~> s-1]
  eps_xy = 0.5*(dudy + dvdx)
  ! dudx and dvdy, which a quarter turn exchanges, are grouped.
  eps_e = sqrt(max(0.0, (((dudx*dudx) + (dvdy*dvdy)) + &
                         ((dudx*dvdy) + (eps_xy*eps_xy)))))
end function dg1_face_eps_eff

!> Effective strain rate at the midpoint of u-face (I,j), from the cell-mean velocities on either
!! side, one-sided where a neighbour is outside the data domain.
pure function dg1_eps_face_u(CS, G, I, j) result(eps_e)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< Ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure
  integer,                intent(in) :: I  !< i-index of the face
  integer,                intent(in) :: j  !< j-index of the face
  real :: eps_e                            !< Effective strain rate [T-1 ~> s-1]
  real :: u_mn, u_pl, v_mn, v_pl ! Cell-mean velocities on the minus/plus sides [L T-1 ~> m s-1]
  real :: dudx, dudy, dvdx, dvdy ! Face-midpoint velocity gradients [T-1 ~> s-1]
  integer :: i_lo, i_hi          ! Neighbour indices clipped to the data domain

  ! The four-corner means are summed in opposite pairs, as every other corner reduction in the
  ! DG path is, because a quarter turn maps each diagonal pair onto the other.  Grouping them by
  ! row, which is the order they are written in, does not survive the turn.
  i_lo = max(I-1, G%isd) ; i_hi = min(I+1, G%ied)
  if (i_lo < I) then
    u_mn = 0.25*((CS%u_shelf(i_lo,j-1) + CS%u_shelf(I,j  )) + &
                 (CS%u_shelf(I   ,j-1) + CS%u_shelf(i_lo,j)))
    v_mn = 0.25*((CS%v_shelf(i_lo,j-1) + CS%v_shelf(I,j  )) + &
                 (CS%v_shelf(I   ,j-1) + CS%v_shelf(i_lo,j)))
  else
    u_mn = 0.5*(CS%u_shelf(I,j-1) + CS%u_shelf(I,j))
    v_mn = 0.5*(CS%v_shelf(I,j-1) + CS%v_shelf(I,j))
  endif
  if (i_hi > I) then
    u_pl = 0.25*((CS%u_shelf(I   ,j-1) + CS%u_shelf(i_hi,j  )) + &
                 (CS%u_shelf(i_hi,j-1) + CS%u_shelf(I   ,j  )))
    v_pl = 0.25*((CS%v_shelf(I   ,j-1) + CS%v_shelf(i_hi,j  )) + &
                 (CS%v_shelf(i_hi,j-1) + CS%v_shelf(I   ,j  )))
  else
    u_pl = 0.5*(CS%u_shelf(I,j-1) + CS%u_shelf(I,j))
    v_pl = 0.5*(CS%v_shelf(I,j-1) + CS%v_shelf(I,j))
  endif
  dudx = (u_pl - u_mn) / G%dxCu(I,j)
  dvdx = (v_pl - v_mn) / G%dxCu(I,j)
  dudy = (CS%u_shelf(I,j) - CS%u_shelf(I,j-1)) / G%dyCu(I,j)
  dvdy = (CS%v_shelf(I,j) - CS%v_shelf(I,j-1)) / G%dyCu(I,j)
  eps_e = dg1_face_eps_eff(dudx, dudy, dvdx, dvdy)
end function dg1_eps_face_u

!> Effective strain rate at the midpoint of v-face (i,J), as dg1_eps_face_u.
pure function dg1_eps_face_v(CS, G, i, J) result(eps_e)
  type(ice_shelf_dyn_CS), intent(in) :: CS !< Ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in) :: G  !< The grid structure
  integer,                intent(in) :: i  !< i-index of the face
  integer,                intent(in) :: J  !< j-index of the face
  real :: eps_e                            !< Effective strain rate [T-1 ~> s-1]
  real :: u_mn, u_pl, v_mn, v_pl ! Cell-mean velocities on the minus/plus sides [L T-1 ~> m s-1]
  real :: dudx, dudy, dvdx, dvdy ! Face-midpoint velocity gradients [T-1 ~> s-1]
  integer :: j_lo, j_hi          ! Neighbour indices clipped to the data domain

  ! Opposite pairs again, matching dg1_eps_face_u so the turn carries one onto the other.
  j_lo = max(J-1, G%jsd) ; j_hi = min(J+1, G%jed)
  if (j_lo < J) then
    u_mn = 0.25*((CS%u_shelf(i-1,j_lo) + CS%u_shelf(i,J   )) + &
                 (CS%u_shelf(i  ,j_lo) + CS%u_shelf(i-1,J  )))
    v_mn = 0.25*((CS%v_shelf(i-1,j_lo) + CS%v_shelf(i,J   )) + &
                 (CS%v_shelf(i  ,j_lo) + CS%v_shelf(i-1,J  )))
  else
    u_mn = 0.5*(CS%u_shelf(i-1,J) + CS%u_shelf(i,J))
    v_mn = 0.5*(CS%v_shelf(i-1,J) + CS%v_shelf(i,J))
  endif
  if (j_hi > J) then
    u_pl = 0.25*((CS%u_shelf(i-1,J   ) + CS%u_shelf(i,j_hi)) + &
                 (CS%u_shelf(i  ,J   ) + CS%u_shelf(i-1,j_hi)))
    v_pl = 0.25*((CS%v_shelf(i-1,J   ) + CS%v_shelf(i,j_hi)) + &
                 (CS%v_shelf(i  ,J   ) + CS%v_shelf(i-1,j_hi)))
  else
    u_pl = 0.5*(CS%u_shelf(i-1,J) + CS%u_shelf(i,J))
    v_pl = 0.5*(CS%v_shelf(i-1,J) + CS%v_shelf(i,J))
  endif
  dudy = (u_pl - u_mn) / G%dyCv(i,J)
  dvdy = (v_pl - v_mn) / G%dyCv(i,J)
  dudx = (CS%u_shelf(i,J) - CS%u_shelf(i-1,J)) / G%dxCv(i,J)
  dvdx = (CS%v_shelf(i,J) - CS%v_shelf(i-1,J)) / G%dxCv(i,J)
  eps_e = dg1_face_eps_eff(dudx, dudy, dvdx, dvdy)
end function dg1_eps_face_v

!> Artificial-viscosity quantities for one face between cells A and B: u_eff and the equivalent
!! thickness jump at each face quadrature point, the face coefficient and jump-mode decay rate,
!! and whether the face is a stagnant jump. Endpoint values are ordered along the face.
subroutine dg1_art_visc_face(CS, bc_A, bc_B, hbc_A_min, hbc_B_min, h_A, h_B, bed, u, v, H_ref, dx_perp, &
                             eps_e_face, rhoi_rhow, advect_grid_inv, advect_inv_L, inv_tau_floor, &
                             ueff, dheq, coef_face, rate_face, idle)
  type(ice_shelf_dyn_CS), intent(in)  :: CS    !< Ice shelf dynamics control structure
  logical,                intent(in)  :: bc_A  !< True if cell A is a thickness boundary (hmask=3)
  logical,                intent(in)  :: bc_B  !< True if cell B is a thickness boundary (hmask=3)
  real,                   intent(in)  :: hbc_A_min !< Lower bound on the face thickness where
                                               !! cell A is a thickness boundary [Z ~> m]
  real,                   intent(in)  :: hbc_B_min !< The same for cell B [Z ~> m]
  real, dimension(2),     intent(in)  :: h_A   !< Cell A corner thicknesses at the face ends [Z ~> m]
  real, dimension(2),     intent(in)  :: h_B   !< Cell B corner thicknesses at the face ends [Z ~> m]
  real, dimension(2),     intent(in)  :: bed   !< Bed at the face ends [Z ~> m]
  real, dimension(2),     intent(in)  :: u     !< Zonal velocity at the face ends [L T-1 ~> m s-1]
  real, dimension(2),     intent(in)  :: v     !< Meridional velocity at the face ends [L T-1 ~> m s-1]
  real,                   intent(in)  :: H_ref !< Reference thickness for the gate [Z ~> m]
  real,                   intent(in)  :: dx_perp !< Across-face length [L ~> m]
  real,                   intent(in)  :: eps_e_face !< Effective strain rate at the face [T-1 ~> s-1]
  real,                   intent(in)  :: rhoi_rhow !< Ice to ocean density ratio [nondim]
  logical,                intent(in)  :: advect_grid_inv !< If true, scale the advective term by dx_perp/L_ref
  real,                   intent(in)  :: advect_inv_L !< 1/DG1_ART_VISC_ADVECT_L_REF, or 0 [L-1 ~> m-1]
  real,                   intent(in)  :: inv_tau_floor !< 1/DG1_ART_VISC_TAU_FLOOR, or 0 [T-1 ~> s-1]
  real, dimension(2),     intent(out) :: ueff  !< u_eff at each face QP [L T-1 ~> m s-1]
  real, dimension(2),     intent(out) :: dheq  !< Equivalent thickness jump at each face QP [Z ~> m]
  real,                   intent(out) :: coef_face !< Face coefficient before the per-cell cap [nondim]
  real,                   intent(out) :: rate_face !< Face jump-mode decay rate [T-1 ~> s-1]
  logical,                intent(out) :: idle  !< True for a stagnant jump

  ! 2-point Gauss-Legendre on [0,1], with 1-gp1 == gp2 and 1-gp2 == gp1 exactly.
  real, parameter :: gp1 = 0.5 * (1.0 - sqrt(1.0/3.0))
  real, parameter :: gp2 = 1.0 - gp1
  real :: t_face, t_co       ! Position along the face and its complement [nondim]
  real :: u_at_qp, v_at_qp   ! Velocity at a face QP [L T-1 ~> m s-1]
  real :: u_mag_qp           ! Speed at a face QP [L T-1 ~> m s-1]
  real :: h_A_qp, h_B_qp     ! Thickness on the two sides at a face QP [Z ~> m]
  real :: bed_qp             ! Bed at a face QP [Z ~> m]
  real :: ds_use             ! Surface jump at a face QP [Z ~> m]
  real :: dh_eq              ! Equivalent thickness jump at a face QP [Z ~> m]
  real :: u_floor_qp         ! Strain-rate term of u_eff [L T-1 ~> m s-1]
  real :: u_adv_qp           ! Advective term of u_eff [L T-1 ~> m s-1]
  real :: u_eff_qp           ! u_eff at a face QP [L T-1 ~> m s-1]
  real :: dh_eq_face_max     ! Max |dh_eq| over the face QPs [Z ~> m]
  real :: ds_face_max        ! Max |ds_use| over the face QPs [Z ~> m]
  real :: u_eff_face_max     ! Max u_eff over the face QPs [L T-1 ~> m s-1]
  real :: u_mag_face_max     ! Max speed over the face QPs [L T-1 ~> m s-1]
  real :: r_face             ! Gate ratio max|[s]|/H_ref [nondim]
  integer :: gp

  dh_eq_face_max = 0.0
  u_mag_face_max = 0.0
  u_eff_face_max = 0.0
  ds_face_max = 0.0
  do gp = 1, 2
    if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
    t_co = 1.0 - t_face
    u_at_qp = (t_co*u(1)) + (t_face*u(2))
    v_at_qp = (t_co*v(1)) + (t_face*v(2))
    u_mag_qp = sqrt((u_at_qp*u_at_qp) + (v_at_qp*v_at_qp))
    ! A boundary side supplies its own nodal face values, exactly as an interior side does, so
    ! that the jump the penalty sees is the real one.  hbc only sets a floor.
    h_A_qp = (t_co*h_A(1)) + (t_face*h_A(2))
    if (bc_A) h_A_qp = max(h_A_qp, hbc_A_min)
    h_B_qp = (t_co*h_B(1)) + (t_face*h_B(2))
    if (bc_B) h_B_qp = max(h_B_qp, hbc_B_min)
    bed_qp = (t_co*bed(1)) + (t_face*bed(2))
    ds_use = dg1_wb_surface_jump(h_A_qp, h_B_qp, bed_qp, rhoi_rhow)
    dh_eq = ds_use * dg1_wb_slope_mean(h_A_qp, h_B_qp, bed_qp, rhoi_rhow)
    u_floor_qp = CS%dg_art_visc_strain_coef * eps_e_face * dx_perp
    if (advect_grid_inv) then
      u_adv_qp = (CS%dg_art_visc_advect_coef * u_mag_qp) * (dx_perp * advect_inv_L)
    else
      u_adv_qp = CS%dg_art_visc_advect_coef * u_mag_qp
    endif
    u_eff_qp = (u_adv_qp + u_floor_qp) + (dx_perp * inv_tau_floor)

    ueff(gp) = u_eff_qp
    dheq(gp) = dh_eq

    dh_eq_face_max = max(dh_eq_face_max, abs(dh_eq))
    ds_face_max = max(ds_face_max, abs(ds_use))
    u_eff_face_max = max(u_eff_face_max, u_eff_qp)
    u_mag_face_max = max(u_mag_face_max, u_mag_qp)
  enddo

  ! Gate on the surface jump relative to the mean thickness.
  r_face = ds_face_max / H_ref
  if (r_face <= 0.0) then
    coef_face = 0.0
  elseif (r_face >= CS%dg_art_visc_r_hi) then
    coef_face = 1.0
  else
    coef_face = r_face / CS%dg_art_visc_r_hi
  endif
  ! Jump-mode decay rate: node localization (x2) and consistent mass (x2) times amp.
  rate_face = 4.0 * DG1_WB_JUMP_RATE_AMP * coef_face * u_eff_face_max / dx_perp

  idle = (u_mag_face_max < CS%dg_slow_idle_u_tiny) .and. &
         (eps_e_face < CS%dg_slow_idle_eps_tiny) .and. &
         (dh_eq_face_max > CS%dg_slow_idle_s_tol)
end subroutine dg1_art_visc_face

!> DG(1) right-hand side M dh/dt = int grad(N).(u h) dA - int N (u.n) h_upwind ds + art visc,
!! with 2x2 Gauss points in the cell and 2 per face. Specified-flux faces split their flux
!! equally between the two corners.
subroutine DG1_nodal_spatial_operator(CS, G, hmask, h_nodal_in, rhs, uh_ice, vh_ice, dt)
  type(ice_shelf_dyn_CS), intent(in) :: CS
  type(ocean_grid_type),  intent(in) :: G
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in)        :: hmask
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(in)    :: h_nodal_in
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(out)   :: rhs
  real, dimension(SZDIB_(G),SZDJ_(G)),    intent(inout) :: uh_ice
  real, dimension(SZDI_(G),SZDJB_(G)),    intent(inout) :: vh_ice
  real,                                   intent(in)    :: dt !< Time step for the art-visc cap [T ~> s]

  ! 2-point Gauss-Legendre on [0,1], with 1-gp1 == gp2 and 1-gp2 == gp1 exactly.
  real, parameter :: gp1 = 0.5 * (1.0 - sqrt(1.0/3.0))
  real, parameter :: gp2 = 1.0 - gp1
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
  real :: visc_flux_qp       ! Artificial-viscosity face flux per QP [Z L2 T-1 ~> m3 s-1]
  real :: Hbar_A, Hbar_B     ! Cell-mean thickness on the two sides of a DG face [Z ~> m]
  real :: H_ref              ! Reference thickness for the gate ratio [Z ~> m]
  real :: coef_face          ! Face viscosity coefficient [nondim]
  real :: scale_AB           ! Smaller cap factor of the two cells [nondim]
  real :: rhoi_rhow_wb       ! Ice/ocean density ratio for the well-balanced jump [nondim]
  real :: advect_inv_L       ! 1/DG1_ART_VISC_ADVECT_L_REF, or 0 [L-1 ~> m-1]
  real :: inv_tau_floor      ! 1/DG1_ART_VISC_TAU_FLOOR, or 0 [T-1 ~> s-1]
  logical :: advect_grid_inv ! If true, scale the advective term by dx_perp/L_ref.
  real :: eps_e_face         ! Effective strain rate at the face midpoint [T-1 ~> s-1]
  logical :: valid_A_visc, valid_B_visc ! Side A/B is hmask 1 or 3.
  logical :: idle_face      ! True for a stagnant-jump face
  ! Art-visc pass-1 results: face coefficient and rate, and per-QP u_eff and Delta h_eq.
  real, dimension(SZDIB_(G),SZDJ_(G))   :: cK_E, rate_E
  real, dimension(SZDI_(G),SZDJB_(G))   :: cK_N, rate_N
  real, dimension(SZDIB_(G),SZDJ_(G),2) :: ueff_E, dheq_E
  real, dimension(SZDI_(G),SZDJB_(G),2) :: ueff_N, dheq_N
  logical, dimension(SZDIB_(G),SZDJ_(G)) :: active_E
  logical, dimension(SZDI_(G),SZDJB_(G)) :: active_N
  real, dimension(SZDI_(G),SZDJ_(G))    :: cell_scale
  real :: S_K                ! Sum of a cell's face rates [T-1 ~> s-1]
  real :: sK_W, sK_E, sK_S, sK_N ! Its four face rates, zero where the face is inactive
                             ! [T-1 ~> s-1]
  real, dimension(2,2,2,2) :: qv_vol ! Per-QP volume contribution to the 4 cell corners,
                                     ! indexed (qx,qy,a,b) [Z L2 T-1 ~> m3 s-1]
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: rhs_vol
  ! Face terms kept per direction and summed pairwise at the end, for rotation invariance.
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: rhs_advx, rhs_advy, rhs_viscx, rhs_viscy

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  rhoi_rhow_wb = CS%density_ice / CS%density_ocean_avg

  rhs(:,:,:,:) = 0.0
  rhs_vol(:,:,:,:) = 0.0
  rhs_advx(:,:,:,:) = 0.0 ; rhs_advy(:,:,:,:) = 0.0
  rhs_viscx(:,:,:,:) = 0.0 ; rhs_viscy(:,:,:,:) = 0.0
  if (associated(CS%dg_art_visc_coef_u)) CS%dg_art_visc_coef_u(:,:) = 0.0
  if (associated(CS%dg_art_visc_coef_v)) CS%dg_art_visc_coef_v(:,:) = 0.0
  if (associated(CS%dg_art_visc_nu_u)) CS%dg_art_visc_nu_u(:,:) = 0.0
  if (associated(CS%dg_art_visc_nu_v)) CS%dg_art_visc_nu_v(:,:) = 0.0
  if (associated(CS%dg_slow_idle_face_u)) CS%dg_slow_idle_face_u(:,:) = 0.0
  if (associated(CS%dg_slow_idle_face_v)) CS%dg_slow_idle_face_v(:,:) = 0.0
  ! 1 means uncapped.
  if (associated(CS%dg_art_visc_cell_scale)) CS%dg_art_visc_cell_scale(:,:) = 1.0

  advect_grid_inv = (CS%dg_art_visc_advect_L_ref > 0.0)
  advect_inv_L = 0.0
  if (advect_grid_inv) advect_inv_L = 1.0 / CS%dg_art_visc_advect_L_ref
  inv_tau_floor = 0.0
  if (CS%dg_art_visc_tau_floor > 0.0) inv_tau_floor = 1.0 / CS%dg_art_visc_tau_floor

  ! Volume integral over each cell.
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    dxCv_S = G%dxCv(i,j-1) ; dxCv_N = G%dxCv(i,j)
    dyCu_W = G%dyCu(i-1,j) ; dyCu_E = G%dyCu(i,j)

    do qy = 1, 2 ; do qx = 1, 2
      if (qx == 1) then ; xi_q  = gp1 ; else ; xi_q  = gp2 ; endif
      if (qy == 1) then ; eta_q = gp1 ; else ; eta_q = gp2 ; endif
      a_qp = (dxCv_S*(1.0 - eta_q)) + (dxCv_N*eta_q)
      d_qp = (dyCu_W*(1.0 - xi_q))  + (dyCu_E*xi_q)

      N11 = (1.0-xi_q)*(1.0-eta_q) ; N21 = xi_q*(1.0-eta_q)
      N12 = (1.0-xi_q)*eta_q       ; N22 = xi_q*eta_q
      dN_dxi_11  = -(1.0 - eta_q) ; dN_dxi_21  =  (1.0 - eta_q)
      dN_dxi_12  = -eta_q         ; dN_dxi_22  =  eta_q
      dN_deta_11 = -(1.0 - xi_q)  ; dN_deta_21 = -xi_q
      dN_deta_12 =  (1.0 - xi_q)  ; dN_deta_22 =  xi_q

      h_qp = (((N11*h_nodal_in(i,j,1,1)) + (N22*h_nodal_in(i,j,2,2))) + &
              ((N21*h_nodal_in(i,j,2,1)) + (N12*h_nodal_in(i,j,1,2))))
      u_qp = (((N11*CS%u_shelf(i-1,j-1)) + (N22*CS%u_shelf(i,j))) + &
              ((N21*CS%u_shelf(i,j-1))   + (N12*CS%u_shelf(i-1,j))))
      v_qp = (((N11*CS%v_shelf(i-1,j-1)) + (N22*CS%v_shelf(i,j))) + &
              ((N21*CS%v_shelf(i,j-1))   + (N12*CS%v_shelf(i-1,j))))

      ! Volume contribution at this qp for each test function N(a,b):
      ! weight * h * (u * dN/dxi * d + v * dN/deta * a), summed below in diagonal QP pairs.
      qv_vol(qx,qy,1,1) = gw*gw * h_qp * (((u_qp * dN_dxi_11) * d_qp) + ((v_qp * dN_deta_11) * a_qp))
      qv_vol(qx,qy,2,1) = gw*gw * h_qp * (((u_qp * dN_dxi_21) * d_qp) + ((v_qp * dN_deta_21) * a_qp))
      qv_vol(qx,qy,1,2) = gw*gw * h_qp * (((u_qp * dN_dxi_12) * d_qp) + ((v_qp * dN_deta_12) * a_qp))
      qv_vol(qx,qy,2,2) = gw*gw * h_qp * (((u_qp * dN_dxi_22) * d_qp) + ((v_qp * dN_deta_22) * a_qp))
    enddo ; enddo

    do b = 1, 2 ; do a = 1, 2
      rhs_vol(i,j,a,b) = (qv_vol(1,1,a,b) + qv_vol(2,2,a,b)) + &
                         (qv_vol(1,2,a,b) + qv_vol(2,1,a,b))
    enddo ; enddo
  enddo ; enddo

  ! East-face fluxes between cells (i,j) and (i+1,j).
  do j = jsc, jec ; do i = isc-1, iec
    if (CS%u_face_mask(i,j) == 4.0) then
      face_flux_total = G%dyCu(i,j) * CS%u_flux_bdry_val(i,j)
      uh_ice(i,j) = uh_ice(i,j) + face_flux_total
      if (i >= isc .and. hmask(i,j) == 1.0) then
        rhs_advx(i,j,2,1) = rhs_advx(i,j,2,1) - (0.5*face_flux_total)
        rhs_advx(i,j,2,2) = rhs_advx(i,j,2,2) - (0.5*face_flux_total)
      endif
      if (i+1 <= iec .and. hmask(i+1,j) == 1.0) then
        rhs_advx(i+1,j,1,1) = rhs_advx(i+1,j,1,1) + (0.5*face_flux_total)
        rhs_advx(i+1,j,1,2) = rhs_advx(i+1,j,1,2) + (0.5*face_flux_total)
      endif
    else if (((i >= isc .and. (hmask(i,j) == 1.0 .or. hmask(i,j) == 3.0))) .or. &
             ((i+1 <= iec .and. (hmask(i+1,j) == 1.0 .or. hmask(i+1,j) == 3.0)))) then
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        u_at_qp = (t_co*CS%u_shelf(i,j-1)) + (t_face*CS%u_shelf(i,j))
        if (u_at_qp >= 0.0) then
          if (hmask(i,j) == 1.0 .or. hmask(i,j) == 3.0) then
            ! A boundary cell carries its own nodal thicknesses, so its face is read the same way
            ! an interior face is.  Reading one number for the whole face would make the boundary
            ! piecewise constant and stop the neighbour holding a gradient against it.
            h_upwind = (t_co*h_nodal_in(i,j,2,1)) + (t_face*h_nodal_in(i,j,2,2))
            if (hmask(i,j) == 3.0) h_upwind = max(h_upwind, CS%min_h_shelf)
          else
            h_upwind = 0.0
          endif
        else
          if (hmask(i+1,j) == 1.0 .or. hmask(i+1,j) == 3.0) then
            h_upwind = (t_co*h_nodal_in(i+1,j,1,1)) + (t_face*h_nodal_in(i+1,j,1,2))
            if (hmask(i+1,j) == 3.0) h_upwind = max(h_upwind, CS%min_h_shelf)
          else
            h_upwind = 0.0
          endif
        endif
        h_upwind = max(h_upwind, 0.0)
        flux_qp = gw * u_at_qp * h_upwind * G%dyCu(i,j)
        uh_ice(i,j) = uh_ice(i,j) + flux_qp
        if (i >= isc .and. hmask(i,j) == 1.0) then
          rhs_advx(i,j,2,1) = rhs_advx(i,j,2,1) - (flux_qp * t_co)
          rhs_advx(i,j,2,2) = rhs_advx(i,j,2,2) - (flux_qp * t_face)
        endif
        if (i+1 <= iec .and. hmask(i+1,j) == 1.0) then
          rhs_advx(i+1,j,1,1) = rhs_advx(i+1,j,1,1) + (flux_qp * t_co)
          rhs_advx(i+1,j,1,2) = rhs_advx(i+1,j,1,2) + (flux_qp * t_face)
        endif
      enddo
    endif
  enddo ; enddo

  ! DG(1) artificial viscosity. Each face QP carries the flux gw*c_face*u_eff*Delta h_eq*ell,
  ! where Delta h_eq is the thickness jump equivalent to the surface jump, so a continuous
  ! surface is not damped. c_face = min(max|[s]|/(H_ref*r_hi), 1). Each cell caps dt times its
  ! summed face rates 4*amp*c*u_eff/dx_perp at DG1_ART_VISC_KCELL. An hmask=3 side is read from
  ! its own corners, like any other side, with min_h_shelf as the only floor.
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

  ! Pass 1, east faces: u_eff and Delta h_eq per QP, c_face and rate per face.
  do j = jsc, jec ; do i = isc-1, iec
    if (CS%u_face_mask(i,j) == 4.0) cycle  ! specified-flux face: viscosity undefined.
    valid_A_visc = (hmask(i,  j) == 1.0 .or. hmask(i,  j) == 3.0)
    valid_B_visc = (hmask(i+1,j) == 1.0 .or. hmask(i+1,j) == 3.0)
    if (.not. (valid_A_visc .and. valid_B_visc)) cycle
    if (.not. (hmask(i,j) == 1.0 .or. hmask(i+1,j) == 1.0)) cycle
    active_E(i,j) = .true.

    if (hmask(i,j) == 3.0) then
      Hbar_A = max(nodal_cell_mean(h_nodal_in(i,j,:,:), CS%cell_mean_w(i,j,:,:)), CS%min_h_shelf)
    else
      Hbar_A = nodal_cell_mean(h_nodal_in(i,  j,:,:), CS%cell_mean_w(i,  j,:,:))
    endif
    if (hmask(i+1,j) == 3.0) then
      Hbar_B = max(nodal_cell_mean(h_nodal_in(i+1,j,:,:), CS%cell_mean_w(i+1,j,:,:)), CS%min_h_shelf)
    else
      Hbar_B = nodal_cell_mean(h_nodal_in(i+1,j,:,:), CS%cell_mean_w(i+1,j,:,:))
    endif
    H_ref = max(CS%min_h_shelf, 0.5*(Hbar_A + Hbar_B))

    eps_e_face = 0.0
    if (CS%dg_art_visc_strain_coef > 0.0) eps_e_face = dg1_eps_face_u(CS, G, i, j)
    call dg1_art_visc_face(CS, hmask(i,j) == 3.0, hmask(i+1,j) == 3.0, &
             CS%min_h_shelf, CS%min_h_shelf, &
             h_nodal_in(i,j,2,:), h_nodal_in(i+1,j,1,:), CS%bed_node(i,j-1:j), &
             CS%u_shelf(i,j-1:j), CS%v_shelf(i,j-1:j), H_ref, G%dxCu(i,j), eps_e_face, &
             rhoi_rhow_wb, advect_grid_inv, advect_inv_L, inv_tau_floor, &
             ueff_E(i,j,:), dheq_E(i,j,:), cK_E(i,j), rate_E(i,j), idle_face)
    if (associated(CS%dg_slow_idle_face_u)) then
      if (idle_face) CS%dg_slow_idle_face_u(i,j) = 1.0
    endif
  enddo ; enddo

  ! North-face fluxes between cells (i,j) and (i,j+1).
  do j = jsc-1, jec ; do i = isc, iec
    if (CS%v_face_mask(i,j) == 4.0) then
      face_flux_total = G%dxCv(i,j) * CS%v_flux_bdry_val(i,j)
      vh_ice(i,j) = vh_ice(i,j) + face_flux_total
      if (j >= jsc .and. hmask(i,j) == 1.0) then
        rhs_advy(i,j,1,2) = rhs_advy(i,j,1,2) - (0.5*face_flux_total)
        rhs_advy(i,j,2,2) = rhs_advy(i,j,2,2) - (0.5*face_flux_total)
      endif
      if (j+1 <= jec .and. hmask(i,j+1) == 1.0) then
        rhs_advy(i,j+1,1,1) = rhs_advy(i,j+1,1,1) + (0.5*face_flux_total)
        rhs_advy(i,j+1,2,1) = rhs_advy(i,j+1,2,1) + (0.5*face_flux_total)
      endif
    else if (((j >= jsc .and. (hmask(i,j) == 1.0 .or. hmask(i,j) == 3.0))) .or. &
             ((j+1 <= jec .and. (hmask(i,j+1) == 1.0 .or. hmask(i,j+1) == 3.0)))) then
      do gp = 1, 2
        if (gp == 1) then ; t_face = gp1 ; else ; t_face = gp2 ; endif
        t_co = 1.0 - t_face
        v_at_qp = (t_co*CS%v_shelf(i-1,j)) + (t_face*CS%v_shelf(i,j))
        if (v_at_qp >= 0.0) then
          if (hmask(i,j) == 1.0 .or. hmask(i,j) == 3.0) then
            h_upwind = (t_co*h_nodal_in(i,j,1,2)) + (t_face*h_nodal_in(i,j,2,2))
            if (hmask(i,j) == 3.0) h_upwind = max(h_upwind, CS%min_h_shelf)
          else
            h_upwind = 0.0
          endif
        else
          if (hmask(i,j+1) == 1.0 .or. hmask(i,j+1) == 3.0) then
            h_upwind = (t_co*h_nodal_in(i,j+1,1,1)) + (t_face*h_nodal_in(i,j+1,2,1))
            if (hmask(i,j+1) == 3.0) h_upwind = max(h_upwind, CS%min_h_shelf)
          else
            h_upwind = 0.0
          endif
        endif
        h_upwind = max(h_upwind, 0.0)
        flux_qp = gw * v_at_qp * h_upwind * G%dxCv(i,j)
        vh_ice(i,j) = vh_ice(i,j) + flux_qp
        if (j >= jsc .and. hmask(i,j) == 1.0) then
          rhs_advy(i,j,1,2) = rhs_advy(i,j,1,2) - (flux_qp * t_co)
          rhs_advy(i,j,2,2) = rhs_advy(i,j,2,2) - (flux_qp * t_face)
        endif
        if (j+1 <= jec .and. hmask(i,j+1) == 1.0) then
          rhs_advy(i,j+1,1,1) = rhs_advy(i,j+1,1,1) + (flux_qp * t_co)
          rhs_advy(i,j+1,2,1) = rhs_advy(i,j+1,2,1) + (flux_qp * t_face)
        endif
      enddo
    endif
  enddo ; enddo

  ! Pass 1, north faces.
  do j = jsc-1, jec ; do i = isc, iec
    if (CS%v_face_mask(i,j) == 4.0) cycle
    valid_A_visc = (hmask(i,j  ) == 1.0 .or. hmask(i,j  ) == 3.0)
    valid_B_visc = (hmask(i,j+1) == 1.0 .or. hmask(i,j+1) == 3.0)
    if (.not. (valid_A_visc .and. valid_B_visc)) cycle
    if (.not. (hmask(i,j) == 1.0 .or. hmask(i,j+1) == 1.0)) cycle
    active_N(i,j) = .true.

    if (hmask(i,j) == 3.0) then
      Hbar_A = max(nodal_cell_mean(h_nodal_in(i,j,:,:), CS%cell_mean_w(i,j,:,:)), CS%min_h_shelf)
    else
      Hbar_A = nodal_cell_mean(h_nodal_in(i,j,  :,:), CS%cell_mean_w(i,j,  :,:))
    endif
    if (hmask(i,j+1) == 3.0) then
      Hbar_B = max(nodal_cell_mean(h_nodal_in(i,j+1,:,:), CS%cell_mean_w(i,j+1,:,:)), CS%min_h_shelf)
    else
      Hbar_B = nodal_cell_mean(h_nodal_in(i,j+1,:,:), CS%cell_mean_w(i,j+1,:,:))
    endif
    H_ref = max(CS%min_h_shelf, 0.5*(Hbar_A + Hbar_B))

    eps_e_face = 0.0
    if (CS%dg_art_visc_strain_coef > 0.0) eps_e_face = dg1_eps_face_v(CS, G, i, j)
    call dg1_art_visc_face(CS, hmask(i,j) == 3.0, hmask(i,j+1) == 3.0, &
             CS%min_h_shelf, CS%min_h_shelf, &
             h_nodal_in(i,j,:,2), h_nodal_in(i,j+1,:,1), CS%bed_node(i-1:i,j), &
             CS%u_shelf(i-1:i,j), CS%v_shelf(i-1:i,j), H_ref, G%dyCv(i,j), eps_e_face, &
             rhoi_rhow_wb, advect_grid_inv, advect_inv_L, inv_tau_floor, &
             ueff_N(i,j,:), dheq_N(i,j,:), cK_N(i,j), rate_N(i,j), idle_face)
    if (associated(CS%dg_slow_idle_face_v)) then
      if (idle_face) CS%dg_slow_idle_face_v(i,j) = 1.0
    endif
  enddo ; enddo

  ! Per-cell cap factor min(1, kcell/(S_K*dt)); non-ice and hmask=3 cells keep 1.
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    ! Opposite faces are summed as a pair, and the two pairs are then added.  A quarter turn
    ! carries the west-east pair onto the south-north pair and back, so it only exchanges the
    ! two pair sums, which addition does not notice.  Summing the four in a fixed E, E, N, N
    ! order does not survive the turn.
    sK_W = 0.0 ; if (active_E(i-1,j)) sK_W = rate_E(i-1,j)
    sK_E = 0.0 ; if (active_E(i  ,j)) sK_E = rate_E(i  ,j)
    sK_S = 0.0 ; if (active_N(i,j-1)) sK_S = rate_N(i,j-1)
    sK_N = 0.0 ; if (active_N(i,j  )) sK_N = rate_N(i,j  )
    S_K = (sK_W + sK_E) + (sK_S + sK_N)
    if (S_K*dt > CS%dg_art_visc_kcell) then
      cell_scale(i,j) = CS%dg_art_visc_kcell / (S_K*dt)
    else
      cell_scale(i,j) = 1.0
    endif
    if (associated(CS%dg_art_visc_cell_scale)) &
      CS%dg_art_visc_cell_scale(i,j) = cell_scale(i,j)
  enddo ; enddo
  ! Both sides of a PE-boundary face need both factors to keep the flux antisymmetric.
  call pass_var(cell_scale, G%domain)

  ! Pass 2, east faces: apply the flux with the smaller of the two cells' cap factors.
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
        rhs_viscx(i,j,2,1) = rhs_viscx(i,j,2,1) + (visc_flux_qp * t_co)
        rhs_viscx(i,j,2,2) = rhs_viscx(i,j,2,2) + (visc_flux_qp * t_face)
      endif
      if (i+1 <= iec .and. hmask(i+1,j) == 1.0) then
        rhs_viscx(i+1,j,1,1) = rhs_viscx(i+1,j,1,1) - (visc_flux_qp * t_co)
        rhs_viscx(i+1,j,1,2) = rhs_viscx(i+1,j,1,2) - (visc_flux_qp * t_face)
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
        rhs_viscy(i,j,1,2) = rhs_viscy(i,j,1,2) + (visc_flux_qp * t_co)
        rhs_viscy(i,j,2,2) = rhs_viscy(i,j,2,2) + (visc_flux_qp * t_face)
      endif
      if (j+1 <= jec .and. hmask(i,j+1) == 1.0) then
        rhs_viscy(i,j+1,1,1) = rhs_viscy(i,j+1,1,1) - (visc_flux_qp * t_co)
        rhs_viscy(i,j+1,2,1) = rhs_viscy(i,j+1,2,1) - (visc_flux_qp * t_face)
      endif
    enddo
  enddo ; enddo

  ! Reduce volume + face.
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    do b = 1, 2 ; do a = 1, 2
      rhs(i,j,a,b) = rhs_vol(i,j,a,b) + &
        ((rhs_advx(i,j,a,b) + rhs_advy(i,j,a,b)) + &
         (rhs_viscx(i,j,a,b) + rhs_viscy(i,j,a,b)))
    enddo ; enddo
  enddo ; enddo
end subroutine DG1_nodal_spatial_operator

!> Biased second difference along one axis: centred (b=1), forward (b=2) or backward (b=3).
!! All give 2*f(0) on the alternating mode and zero on a uniform field.
pure function dg_d2_biased(f, b) result(A)
  real, dimension(-4:4), intent(in) :: f !< Samples along the axis [Z ~> m]
  integer,               intent(in) :: b !< 1 centred, 2 forward, 3 backward
  real :: A                              !< Second difference [Z ~> m]
  select case (b)
    case (1) ; A = f(0) - (0.5*(f(-1) + f(1)))
    case (2) ; A = 0.5*(f(0) - (2.0*f(1))+ f(2))
    case default ; A = 0.5*(f(0) - (2.0*f(-1))+ f(-2))
  end select
end function dg_d2_biased

!> Biased first-difference stencil weights on offsets -2..2, for one axis.
pure function dg_d1_weights(b) result(w)
  integer, intent(in) :: b    !< Bias along the axis: 1 centred, 2 forward, 3 backward
  real, dimension(-2:2) :: w  !< First-difference stencil weights [nondim]
  ! Second order and centred on the cell, so the reference sits where the detector does.
  w(:) = 0.0
  select case (b)
    case (1) ; w(-1) = -0.5 ; w(1) = 0.5
    case (2) ; w(0) = -1.5 ; w(1) = 2.0 ; w(2) = -0.5
    case default ; w(0) = 1.5 ; w(-1) = -2.0 ; w(-2) = 0.5
  end select
end function dg_d1_weights

!> Mean-supported twist at offset (p,q): the cross difference of the cell means, each axis
!! with its own bias. It plays the role of A_ref for the twist.
pure function dg_wref_at(hw, p, q, bx, by) result(wr)
  real, dimension(-4:4,-4:4), intent(in) :: hw !< Cell means on the gather [Z ~> m]
  integer, intent(in) :: p  !< Offset along xi
  integer, intent(in) :: q  !< Offset along eta
  integer, intent(in) :: bx !< Bias along xi
  integer, intent(in) :: by !< Bias along eta
  real :: wr                !< Mean-supported twist [Z ~> m]
  real, dimension(-2:2) :: wxi, weta   ! Per-axis first-difference weights [nondim]
  real, dimension(-2:2,-2:2) :: tm     ! Per-offset contribution to the twist [Z ~> m]
  integer :: m, n, r
  ! One offset from each of the six four-member orbits of (m,n) -> (n,-m) on the 5x5 stencil.
  integer, dimension(6), parameter :: rep_m = (/ 1, 1, 2, 2, 1, 2 /)
  integer, dimension(6), parameter :: rep_n = (/ 0, 1, 0, 1, 2, 2 /)

  ! One weighted stencil summed orbit by orbit in opposite pairs, for rotation invariance.
  wxi = dg_d1_weights(bx) ; weta = dg_d1_weights(by)
  tm(:,:) = 0.0
  do n = -2, 2
    if (weta(n) == 0.0) cycle
    do m = -2, 2
      if (wxi(m) == 0.0) cycle
      tm(m,n) = (wxi(m)*weta(n)) * hw(p+m, q+n)
    enddo
  enddo
  wr = tm(0,0)
  do r = 1, 6
    m = rep_m(r) ; n = rep_n(r)
    wr = wr + ((tm(m,n) + tm(-m,-n)) + (tm(n,-m) + tm(-n,m)))
  enddo
end function dg_wref_at

!> One-dimensional tilt detector A (a second difference of the cell tilt), with the same operator
!! on the bed (A_bed) and on the tilt implied by the cell means (A_ref). The stencil is centred
!! where possible and biased forward or backward next to unusable cells; all forms give 2*t on
!! the alternating mode. The reference is second order at the cell centre at every offset, so
!! A - A_ref vanishes faster than A on a smooth solution. This needs a +-4 gather.
pure subroutine dg_tilt_detector_1d(tv, bv, hv, ok, fv, kink, fit_tol, A, A_bed, A_ref, href, &
                                    valid, mode, kink_c)
  real,    dimension(-4:4), intent(in)  :: tv !< Per-cell tilt along the stencil [Z ~> m]
  real,    dimension(-4:4), intent(in)  :: bv !< Per-cell bed tilt along the stencil [Z ~> m]
  real,    dimension(-4:4), intent(in)  :: hv !< Per-cell mean thickness along the stencil [Z ~> m]
  logical, dimension(-4:4), intent(in)  :: ok !< True where the cell is usable
  real,    dimension(-4:4), intent(in)  :: fv !< Grounded fraction on the gather [nondim]
  logical, intent(in) :: kink !< If true, build the reference with the flotation break in it
  real,    intent(in) :: fit_tol !< Relative slope misfit at which the two-branch model is
                                !! disbelieved, or non-positive to skip the test [nondim]
  real,    intent(out) :: A     !< Tilt-Laplacian detector [Z ~> m]
  real,    intent(out) :: A_bed !< The same operator on the bed [Z ~> m]
  real,    intent(out) :: A_ref !< The same operator on the mean-supported tilt [Z ~> m]
  real,    intent(out) :: href  !< Mean thickness over the stencil used [Z ~> m]
  logical, intent(out) :: valid !< False if no admissible stencil exists
  integer, intent(out) :: mode  !< Stencil chosen: 1 centred, 2 forward, 3 backward, 0 none
  real,    intent(out) :: kink_c  !< Confidence in the kink in the reference, 0 to 1 [nondim]

  real :: rm, r0, rp  ! Mean-supported tilt at the three stencil cells [Z ~> m]
  real :: s_g, s_f    ! Branch slopes across a flotation transition [Z ~> m]
  real :: sigma       ! How much of a transition the reference's window spans [nondim]
  real :: cf          ! Confidence in the branch slopes that were found [nondim]
  integer :: km, k0, kp ! Offsets of the three cells whose tilt the reference supplies
  integer :: rlo, rhi   ! Offsets of the cell means the reference reads
  logical :: kv       ! True if the branch slopes could be built

  A = 0.0 ; A_bed = 0.0 ; A_ref = 0.0 ; href = 0.0 ; valid = .true. ; mode = 1
  kink_c = 0.0

  if (all(ok(-2:2))) then                       ! centred
    A     = tv(0) - (0.5*(tv(-1) + tv(1)))
    A_bed = bv(0) - (0.5*(bv(-1) + bv(1)))
    rm = 0.5*(hv(0) - hv(-2)) ; r0 = 0.5*(hv(1) - hv(-1)) ; rp = 0.5*(hv(2) - hv(0))
    km = -1 ; k0 = 0 ; kp = 1 ; rlo = -2 ; rhi = 2
    href  = (hv(0) + (hv(-1) + hv(1))) / 3.0
  elseif (all(ok(0:4))) then                    ! forward
    A     = 0.5*(tv(0) - (2.0*tv(1))+ tv(2))
    A_bed = 0.5*(bv(0) - (2.0*bv(1))+ bv(2))
    r0 = 0.5*((-3.0*hv(0))+ (4.0*hv(1))- hv(2))
    rm = 0.5*((-3.0*hv(1))+ (4.0*hv(2))- hv(3))
    rp = 0.5*((-3.0*hv(2))+ (4.0*hv(3))- hv(4))
    km = 1 ; k0 = 0 ; kp = 2 ; rlo = 0 ; rhi = 4
    href  = ((hv(0) + hv(1)) + hv(2)) / 3.0 ; mode = 2
  elseif (all(ok(-4:0))) then                   ! backward
    A     = 0.5*(tv(0) - (2.0*tv(-1))+ tv(-2))
    A_bed = 0.5*(bv(0) - (2.0*bv(-1))+ bv(-2))
    r0 = 0.5*( (3.0*hv(0))- (4.0*hv(-1))+ hv(-2))
    rm = 0.5*( (3.0*hv(-1))- (4.0*hv(-2))+ hv(-3))
    rp = 0.5*( (3.0*hv(-2))- (4.0*hv(-3))+ hv(-4))
    km = -1 ; k0 = 0 ; kp = -2 ; rlo = -4 ; rhi = 0
    href  = ((hv(0) + hv(-1)) + hv(-2)) / 3.0 ; mode = 3
  else
    valid = .false. ; mode = 0
  endif
  if (.not.valid) return

  ! Blend in the flotation kink by the span of the grounded fraction; zero span changes nothing.
  if (kink) then
    call dg_kink_branches(hv, fv, ok, rlo, rhi, fit_tol, s_g, s_f, sigma, cf, kv)
    if (kv) then
      sigma = sigma * cf
      rm = ((1.0-sigma)*rm) + (sigma*((fv(km)*s_g)+ ((1.0-fv(km))*s_f)))
      r0 = ((1.0-sigma)*r0) + (sigma*((fv(k0)*s_g)+ ((1.0-fv(k0))*s_f)))
      rp = ((1.0-sigma)*rp) + (sigma*((fv(kp)*s_g)+ ((1.0-fv(kp))*s_f)))
      kink_c = cf
    endif
  endif

  if (mode == 1) then
    A_ref = r0 - (0.5*(rm + rp))
  else
    A_ref = 0.5*(r0 - (2.0*rm)+ rp)
  endif
end subroutine dg_tilt_detector_1d

!> Grounded and floating branch slopes across a flotation transition, from the cells nearest it
!! lying wholly on each side, and the span of the grounded fraction. A piecewise-linear profile
!! breaking at fraction f of a cell rises f*s_g + (1-f)*s_f across it. valid is false if either
!! side has fewer than two such cells.
pure subroutine dg_kink_branches(hv, fv, ok, lo, hi, fit_tol, s_g, s_f, sigma, conf, valid)
  real,    dimension(-4:4), intent(in)  :: hv !< Cell means on the gather [Z ~> m]
  real,    dimension(-4:4), intent(in)  :: fv !< Grounded fraction on the gather [nondim]
  logical, dimension(-4:4), intent(in)  :: ok !< True where the cell is usable
  integer, intent(in)  :: lo    !< First offset the reference reads
  integer, intent(in)  :: hi    !< Last offset the reference reads
  real,    intent(in)  :: fit_tol !< Relative slope misfit giving zero confidence, or
                                !! non-positive to skip the fit test [nondim]
  real,    intent(out) :: s_g   !< Rise per cell on the grounded branch [Z ~> m]
  real,    intent(out) :: s_f   !< Rise per cell on the floating branch [Z ~> m]
  real,    intent(out) :: sigma !< Span max(f) - min(f) over the window [nondim]
  real,    intent(out) :: conf  !< How far the branch slopes are to be trusted, 0 to 1 [nondim]
  logical, intent(out) :: valid !< False if either branch has fewer than two cells
  real, parameter :: eps = 1.0e-9 ! Rounding tolerance on wholly grounded or afloat [nondim]
  integer :: k, ng, nf, ge, fs, slo, shi
  logical :: gnd_low    ! True when the grounded end of the window is the low one
  real :: den           ! Largest slope in play, the scale a misfit is judged against [Z ~> m]
  real :: amax          ! Worst centre-to-centre misfit of the two-branch model [Z ~> m]
  real :: mk, mk1       ! Model rise across cells k and k+1 [Z ~> m]

  s_g = 0.0 ; s_f = 0.0 ; sigma = 0.0 ; conf = 0.0 ; valid = .false.
  sigma = maxval(fv(lo:hi)) - minval(fv(lo:hi))
  if (.not.all(ok(lo:hi))) return

  ! Branch cells are sought over the whole usable run containing the host.
  slo = 0 ; do k = 0, -4, -1 ; if (ok(k)) then ; slo = k ; else ; exit ; endif ; enddo
  shi = 0 ; do k = 0, 4       ; if (ok(k)) then ; shi = k ; else ; exit ; endif ; enddo
  gnd_low = (fv(slo) >= fv(shi))

  ng = 0 ; nf = 0
  if (gnd_low) then
    do k = slo, shi ; if (fv(k) >= 1.0-eps) then ; ng = ng + 1 ; else ; exit ; endif ; enddo
    do k = shi, slo, -1 ; if (fv(k) <= eps) then ; nf = nf + 1 ; else ; exit ; endif ; enddo
    ge = slo + ng - 1 ; fs = shi - nf + 1
  else
    do k = shi, slo, -1 ; if (fv(k) >= 1.0-eps) then ; ng = ng + 1 ; else ; exit ; endif ; enddo
    do k = slo, shi ; if (fv(k) <= eps) then ; nf = nf + 1 ; else ; exit ; endif ; enddo
    ge = shi - ng + 1 ; fs = slo + nf - 1
  endif
  if ((ng < 2) .or. (nf < 2)) return
  ! A two-cell branch uses a one-sided difference, which can pick up an alternating mean.
  conf = 1.0
  if (min(ng, nf) < 3) conf = 0.7

  ! Slopes in the +index sense, nearest the transition.
  if (gnd_low) then
    if (ng >= 3) then ; s_g = 0.5*(hv(ge) - hv(ge-2)) ; else ; s_g = hv(ge) - hv(ge-1) ; endif
    if (nf >= 3) then ; s_f = 0.5*(hv(fs+2) - hv(fs)) ; else ; s_f = hv(fs+1) - hv(fs) ; endif
  else
    if (ng >= 3) then ; s_g = 0.5*(hv(ge+2) - hv(ge)) ; else ; s_g = hv(ge+1) - hv(ge) ; endif
    if (nf >= 3) then ; s_f = 0.5*(hv(fs) - hv(fs-2)) ; else ; s_f = hv(fs) - hv(fs-1) ; endif
  endif

  ! Fit test: compare the model's centre-to-centre rises with the cell means, relative to the
  ! largest slope. Confidence can only fall.
  if (fit_tol > 0.0) then
    den = max(abs(s_g - s_f), max(abs(s_g), abs(s_f)))
    if (den > 0.0) then
      amax = 0.0
      do k = lo, hi-1
        mk  = (fv(k)  *s_g) + ((1.0 - fv(k))  *s_f)
        mk1 = (fv(k+1)*s_g) + ((1.0 - fv(k+1))*s_f)
        amax = max(amax, abs((hv(k+1) - hv(k)) - (0.5*(mk + mk1))))
      enddo
      conf = min(conf, max(0.0, 1.0 - ((amax/den)/fit_tol)))
    endif
  endif
  valid = .true.
end subroutine dg_kink_branches

!> Exact area fraction of the unit cell where a linear function, given at the corners, is positive.
pure function dg_lin_area(d) result(a)
  real, dimension(4), intent(in) :: d !< Corner values, SW SE NW NE [Z ~> m]
  real :: a                           !< Area fraction where the value is positive [nondim]
  real, dimension(2,8) :: p   ! Clipped polygon vertices, cell coordinates [nondim]
  real, dimension(2,4) :: v   ! The square, counter-clockwise [nondim]
  real, dimension(4)   :: dv  ! The function at those vertices, in the same order [Z ~> m]
  real :: t
  integer :: k, kn, n
  v(:,1) = (/ 0.0, 0.0 /) ; v(:,2) = (/ 1.0, 0.0 /)
  v(:,3) = (/ 1.0, 1.0 /) ; v(:,4) = (/ 0.0, 1.0 /)
  dv = (/ d(1), d(2), d(4), d(3) /)     ! SW SE NE NW, matching v
  n = 0
  do k = 1, 4
    kn = 1 + mod(k, 4)
    if (dv(k) > 0.0) then
      n = n + 1 ; p(:,n) = v(:,k)
    endif
    if ((dv(k) > 0.0) .neqv. (dv(kn) > 0.0)) then
      t = dv(k) / (dv(k) - dv(kn))
      n = n + 1 ; p(:,n) = v(:,k) + (t*(v(:,kn) - v(:,k)))
    endif
  enddo
  a = 0.0
  do k = 1, n
    kn = 1 + mod(k, n)
    a = a + ((p(1,k)*p(2,kn))- (p(1,kn)*p(2,k)))
  enddo
  a = min(max(0.5*abs(a), 0.0), 1.0)
end function dg_lin_area

!> Constant to add to a linear function so its positive area fraction is f, by 50 bisection steps.
pure function dg_shift_to_frac(d, f) result(c)
  real, dimension(4), intent(in) :: d !< Corner values before the shift [Z ~> m]
  real, intent(in) :: f               !< Target positive-area fraction [nondim]
  real :: c                           !< Constant to add to every corner [Z ~> m]
  real :: lo, hi, mid
  integer :: it
  if (f <= 0.0) then ; c = -maxval(d) - 1.0 ; return ; endif
  if (f >= 1.0) then ; c = -minval(d) + 1.0 ; return ; endif
  lo = -maxval(d) ; hi = -minval(d)
  do it = 1, 50
    mid = 0.5*(lo + hi)
    if (dg_lin_area(d + mid) < f) then ; lo = mid ; else ; hi = mid ; endif
  enddo
  c = 0.5*(lo + hi)
end function dg_shift_to_frac

!> The flotation break near a cell, as one clipped linear function.
!!
!! The twist needs more than the tilt did.  The tilt depends only on the grounded
!! AREA fraction, so one number placed it; the twist depends on the ORIENTATION
!! as well.  Continuity along a break forces the gradient jump to be normal to it,
!! so a break aligned with x cannot change the x-slope and produces exactly zero
!! twist, while a diagonal one produces the most -- at the same area fraction.
!!
!! Writing the kink as a single clipped linear function carries the orientation
!! for free.  The field near the break is P_f + max(0, D) with D = P_g - P_f
!! linear and vanishing on the break, so D's gradient IS the jump and its zero
!! contour IS the grounding line.  Nothing needs integrating: a cell's tilt and
!! twist follow from D at its four corners.  The 1D rule of dg_kink_branches is
!! this one's shadow, so the two agree by construction.
!!
!! Only the GRADIENT of D has to be estimated, from cells lying wholly on each
!! side; its constant comes from requiring the grounded area of one straddling
!! cell to match ground_frac, which keeps the reference agreeing with the friction
!! and the driving stress about where the grounding line is.  One line is placed
!! for the whole neighbourhood rather than one per cell, a grounding line being
!! one curve and not a set of independent segments.
pure subroutine dg_kink_plane_2d(hw, fw, okw, lo, hi, tol, dqx, dqy, dq0, sigma, &
                                 conf, valid)
  real,    dimension(-4:4,-4:4), intent(in) :: hw  !< Cell means on the gather [Z ~> m]
  real,    dimension(-4:4,-4:4), intent(in) :: fw  !< Grounded fraction [nondim]
  logical, dimension(-4:4,-4:4), intent(in) :: okw !< True where the cell is usable
  integer, intent(in)  :: lo, hi !< Offsets, both axes, over which the span is measured
  real,    intent(out) :: dqx    !< d(D)/dx, one cell [Z ~> m]
  real,    intent(out) :: dqy    !< d(D)/dy, one cell [Z ~> m]
  real,    intent(out) :: dq0    !< Constant, so that D = dq0 + dqx*x + dqy*y [Z ~> m]
  real,    intent(out) :: sigma  !< Span of the grounded fraction over the window [nondim]
  real,    intent(in)  :: tol    !< Misfit at which confidence in the fit reaches zero [nondim]
  real,    intent(out) :: conf   !< How far the single-line fit is to be trusted, 0 to 1 [nondim]
  logical, intent(out) :: valid  !< False if either branch or the anchor is missing
  real, parameter :: eps = 1.0e-9
  real :: gx(2), gy(2), best, dr(4), mis, amax
  integer :: p, q, b, nx(2), ny(2), pa, qa
  logical :: m1, m2

  dqx = 0.0 ; dqy = 0.0 ; dq0 = 0.0 ; conf = 0.0 ; valid = .false.
  sigma = maxval(fw(lo:hi,lo:hi)) - minval(fw(lo:hi,lo:hi))

  ! Branch gradients averaged over all wholly grounded or floating adjacent pairs.
  gx = 0.0 ; gy = 0.0 ; nx = 0 ; ny = 0
  do b = 1, 2
    do q = -4, 4 ; do p = -4, 3
      m1 = okw(p,q) .and. okw(p+1,q)
      if (b == 1) then
        m1 = m1 .and. (fw(p,q) >= 1.0-eps) .and. (fw(p+1,q) >= 1.0-eps)
      else
        m1 = m1 .and. (fw(p,q) <= eps) .and. (fw(p+1,q) <= eps)
      endif
      if (m1) then ; gx(b) = gx(b) + (hw(p+1,q) - hw(p,q)) ; nx(b) = nx(b) + 1 ; endif
    enddo ; enddo
    do q = -4, 3 ; do p = -4, 4
      m2 = okw(p,q) .and. okw(p,q+1)
      if (b == 1) then
        m2 = m2 .and. (fw(p,q) >= 1.0-eps) .and. (fw(p,q+1) >= 1.0-eps)
      else
        m2 = m2 .and. (fw(p,q) <= eps) .and. (fw(p,q+1) <= eps)
      endif
      if (m2) then ; gy(b) = gy(b) + (hw(p,q+1) - hw(p,q)) ; ny(b) = ny(b) + 1 ; endif
    enddo ; enddo
  enddo
  if (any(nx == 0) .or. any(ny == 0)) return
  dqx = (gx(1)/real(nx(1))) - (gx(2)/real(nx(2)))
  dqy = (gy(1)/real(ny(1))) - (gy(2)/real(ny(2)))
  if ((abs(dqx) + abs(dqy)) <= 0.0) return

  ! Anchor the constant on the most straddled cell.
  best = -1.0 ; pa = 0 ; qa = 0
  do q = lo, hi ; do p = lo, hi
    if (.not.okw(p,q)) cycle
    if (0.5 - abs(fw(p,q) - 0.5) > best) then
      best = 0.5 - abs(fw(p,q) - 0.5) ; pa = p ; qa = q
    endif
  enddo ; enddo
  if (best <= 0.0) return

  dr(1) = (dqx*real(pa))       + (dqy*real(qa))
  dr(2) = (dqx*real(pa+1))     + (dqy*real(qa))
  dr(3) = (dqx*real(pa))       + (dqy*real(qa+1))
  dr(4) = (dqx*real(pa+1))     + (dqy*real(qa+1))
  dq0 = dg_shift_to_frac(dr, fw(pa,qa))

  ! Does one line actually describe this neighbourhood?  Having placed it,
  ! predict every straddling cell's grounded fraction from it and compare with
  ! what the partition reported.  A single line that fits reproduces them all.
  ! This is not a formality.  Against a grounding line curving with a radius of
  ! four to eight cells the single-line reference is WORSE than the
  ! mean-supported one it replaces -- errors above the amplitude of the mode it
  ! is meant to isolate, where the old reference's error is roughly constant
  ! because it makes no geometric assumption to be wrong about.  The misfit
  ! separated the two outcomes cleanly in that test, at most 0.16 wherever the
  ! kink model won and at least 0.23 wherever it lost, so confidence is ramped
  ! down over it rather than switched, and at zero confidence the reference is
  ! exactly the one that was there before.
  amax = 0.0
  do q = lo, hi ; do p = lo, hi
    if (.not.okw(p,q)) cycle
    if ((fw(p,q) <= eps) .or. (fw(p,q) >= 1.0-eps)) cycle
    dr(1) = dq0 + (dqx*real(p))   + (dqy*real(q))
    dr(2) = dq0 + (dqx*real(p+1)) + (dqy*real(q))
    dr(3) = dq0 + (dqx*real(p))   + (dqy*real(q+1))
    dr(4) = dq0 + (dqx*real(p+1)) + (dqy*real(q+1))
    mis = abs(dg_lin_area(dr) - fw(p,q))
    amax = max(amax, mis)
  enddo ; enddo
  conf = max(0.0, 1.0 - (amax / max(tol, tiny(1.0))))
  ! Valid means the branches were found, not that the fit is any good: a fit
  ! that is credible but poor must still hand back a confidence for the caller
  ! to interpolate on, rather than being reported as absent.
  valid = .true.
end subroutine dg_kink_plane_2d

!> Twist of max(0, D) over cell (p,q); the linear part has none.
!! The smooth part of the field is a plane and contributes no twist at all, so the
!! whole of it comes from the clipped part and follows from four corner values.
pure function dg_kink_twist_at(dqx, dqy, dq0, p, q) result(w)
  real,    intent(in) :: dqx, dqy, dq0 !< The linear function D [Z ~> m]
  integer, intent(in) :: p, q          !< Cell offsets
  real :: w                            !< Twist the kink implies [Z ~> m]
  real :: m(4)
  m(1) = max(0.0, dq0 + (dqx*real(p))   + (dqy*real(q)))
  m(2) = max(0.0, dq0 + (dqx*real(p+1)) + (dqy*real(q)))
  m(3) = max(0.0, dq0 + (dqx*real(p))   + (dqy*real(q+1)))
  m(4) = max(0.0, dq0 + (dqx*real(p+1)) + (dqy*real(q+1)))
  w = (m(4) - m(3)) - (m(2) - m(1))
end function dg_kink_twist_at

!> The part of detector A outside the interval spanned by 0 and the references r1 and r2.
!! It has the sign of A and never exceeds it, and equals A when both references are zero.
pure function dg_unexplained(A, r1, r2) result(Ad)
  real, intent(in) :: A   !< The detector reading [Z ~> m]
  real, intent(in) :: r1  !< The same operator on the mean-supported tilt [Z ~> m]
  real, intent(in) :: r2  !< The same operator on the bed-forced tilt [Z ~> m]
  real :: Ad              !< The part of A no reference accounts for [Z ~> m]
  real :: lo, hi          ! Ends of the interval the references span [Z ~> m]
  lo = min(min(r1, r2), 0.0) ; hi = max(max(r1, r2), 0.0)
  Ad = A - min(max(A, lo), hi)
end function dg_unexplained

!> Grounding-line protection weight: 0 if the cell (reach 0), the detector stencil (reach 1) or
!! the reference stencil (reach 2) contains both grounded and floating ice, otherwise 1.
pure function dg_gl_reach_wt(fv, mode, reach) result(wt)
  real, dimension(-4:4), intent(in) :: fv !< Grounded fraction on the gather [nondim]
  integer, intent(in) :: mode  !< Stencil chosen: 1 centred, 2 forward, 3 backward
  integer, intent(in) :: reach !< 0 the cell alone, 1 the reach of A, 2 the reach of A_ref
  real :: wt                   !< Weight to apply to the rate [nondim]
  ! Offsets read by A and A_ref for each stencil bias.
  integer, dimension(3), parameter :: Alo = (/ -1, 0, -2 /), Ahi = (/ 1, 2, 0 /)
  integer, dimension(3), parameter :: Rlo = (/ -2, 0, -4 /), Rhi = (/ 2, 4, 0 /)
  integer :: lo, hi
  if (reach <= 0) then
    wt = merge(0.0, 1.0, (fv(0) > 0.0) .and. (fv(0) < 1.0))
    return
  elseif (reach == 1) then
    lo = Alo(mode) ; hi = Ahi(mode)
  else
    lo = Rlo(mode) ; hi = Rhi(mode)
  endif
  wt = merge(0.0, 1.0, (minval(fv(lo:hi)) < 1.0) .and. (maxval(fv(lo:hi)) > 0.0))
end function dg_gl_reach_wt

!> Fraction of full strength that the mode damper runs at, from one of three gates. The
!! thickness gate asks how large the zigzag is against the ice thickness; it closes on thick
!! ice, and because the detector output scales with the cell size while the thickness does not,
!! the slope error it accepts doubles each time the grid is refined by two. The slope gate asks
!! how steep the zigzag is against the real surface slope; both scale with the cell size, so it
!! carries no grid spacing, but it closes on steep grounded ice. The agreement gate asks how
!! much the detector readings agree, which carries no thickness, no slope and no grid spacing.
pure function dg_damp_gate_frac(form, dsdh, excess, spread, href, dsurf, r_hi, rho_s, &
                                floor_rise, rho_g) result(gam)
  integer, intent(in) :: form   !< Gate form, one of DAMP_GATE_THICKNESS, DAMP_GATE_SLOPE
                                !! or DAMP_GATE_AGREEMENT
  real,    intent(in) :: dsdh   !< Surface elevation change per thickness change [nondim]
  real,    intent(in) :: excess !< Detector output the damper would act on [Z ~> m]
  real,    intent(in) :: spread !< Largest minus smallest detector reading [Z ~> m]
  real,    intent(in) :: href   !< Reference thickness for the thickness gate [Z ~> m]
  real,    intent(in) :: dsurf  !< Real rise of the surface across the cell [Z ~> m]
  real,    intent(in) :: r_hi   !< Zigzag over thickness at which the thickness gate opens [nondim]
  real,    intent(in) :: rho_s  !< Zigzag slope over real slope at which the slope gate opens [nondim]
  real,    intent(in) :: floor_rise !< Smallest surface rise the slope gate uses [Z ~> m]
  real,    intent(in) :: rho_g  !< Disagreement over agreement at which the agreement gate is
                                !! half open [nondim]
  real :: gam                   !< Gate fraction, 0 to 1 [nondim]

  real :: denom  ! Reference the detector output is compared with [Z ~> m]

  gam = 0.0
  if (excess <= 0.0) return

  if (form == DAMP_GATE_AGREEMENT) then
    ! Scale-free in amplitude: a pure zigzag gives equal readings, hence no spread and a gate
    ! of 1, whatever the thickness, the slope or the size of the zigzag.
    gam = excess / (excess + (rho_g * max(spread, 0.0)))
  else
    if (form == DAMP_GATE_SLOPE) then
      denom = rho_s * max(abs(dsurf), floor_rise)
    else
      denom = r_hi * href
    endif
    if (denom <= 0.0) return
    gam = min(1.0, (dsdh*excess) / denom)
  endif
end function dg_damp_gate_frac

!> Part of a set of readings that every member agrees on, and the width of the set. The minmod
!! is zero unless all the readings share a sign, and is otherwise the one of smallest magnitude,
!! so the damper never removes more than the least any reading supports. Both outputs come from
!! minval and maxval over the whole set, which do not depend on its order, so a rotation or a
!! reflection of the grid cannot change them.
pure subroutine dg_minmod_spread(nrd, rd, A, S)
  integer,              intent(in)  :: nrd !< How many readings are set
  real, dimension(nrd), intent(in)  :: rd  !< The readings [Z ~> m]
  real,                 intent(out) :: A   !< The part they all agree on [Z ~> m]
  real,                 intent(out) :: S   !< Largest minus smallest reading [Z ~> m]

  real :: rmin, rmax ! Extremes of the set [Z ~> m]
  integer :: n

  A = 0.0 ; S = 0.0
  if (nrd < 1) return
  ! One pass rather than two.  The extremes of a set do not depend on the order in which a
  ! quarter turn presents it, so this reduction needs no pairing.
  rmin = rd(1) ; rmax = rd(1)
  do n = 2, nrd
    rmin = min(rmin, rd(n)) ; rmax = max(rmax, rd(n))
  enddo
  S = rmax - rmin
  if (rmin > 0.0) then
    A = rmin
  elseif (rmax < 0.0) then
    A = rmax
  endif
end subroutine dg_minmod_spread

!> Readings of the three 3-cell stencils that contain a cell, on one field, appended to a set.
!! A reading is the second difference of the field per unit length, scaled back by the width of
!! the centre cell, which is what gives zero for a constant slope on cells of unequal width. The
!! left and right forms are the same expression with mirrored indices, and the centred form pairs
!! the two neighbours, so a rotation or a reflection maps one reading onto another exactly.
!> Weights that turn the cell means of a short run into the tilt those means imply at one cell
!! of it.  The run holds n cells of equal width, numbered 1 to n, and the host is cell p+1.
!! The polynomial of degree n-1 whose cell means match is taken, and its tilt at the host is
!! returned as w . hbar.  The values are exact rationals, so the compiler folds them.
pure function dg_ref_w(n, p) result(w)
  integer, intent(in) :: n   !< Cells in the run, 2 to 5
  integer, intent(in) :: p   !< Position of the host in the run, 0 to n-1
  real, dimension(5) :: w    !< The weights, the last 5-n of them zero [nondim]

  w(:) = 0.0
  select case (10*n + p)
    case (20, 21)
      w(1:2) = (/ -1.0, 1.0 /)
    case (30)
      w(1:3) = (/ -1.5, 2.0, -0.5 /)
    case (31)
      w(1:3) = (/ -0.5, 0.0, 0.5 /)
    case (32)
      w(1:3) = (/ 0.5, -2.0, 1.5 /)
    case (40)
      w(1:4) = (/ -109.0/60.0, 59.0/20.0, -29.0/20.0, 19.0/60.0 /)
    case (41)
      w(1:4) = (/ -19.0/60.0, -11.0/20.0, 21.0/20.0, -11.0/60.0 /)
    case (42)
      w(1:4) = (/ 11.0/60.0, -21.0/20.0, 11.0/20.0, 19.0/60.0 /)
    case (43)
      w(1:4) = (/ -19.0/60.0, 29.0/20.0, -59.0/20.0, 109.0/60.0 /)
    case (50)
      w(1:5) = (/ -49.0/24.0, 77.0/20.0, -14.0/5.0, 73.0/60.0, -9.0/40.0 /)
    case (51)
      w(1:5) = (/ -9.0/40.0, -11.0/12.0, 8.0/5.0, -11.0/20.0, 11.0/120.0 /)
    case (52)
      w(1:5) = (/ 11.0/120.0, -41.0/60.0, 0.0, 41.0/60.0, -11.0/120.0 /)
    case (53)
      w(1:5) = (/ -11.0/120.0, 11.0/20.0, -8.0/5.0, 11.0/12.0, 9.0/40.0 /)
    case (54)
      w(1:5) = (/ 9.0/40.0, -73.0/60.0, 14.0/5.0, -77.0/20.0, 49.0/24.0 /)
  end select
end function dg_ref_w

!> Estimates of the mode at a cell that has lost a neighbour, from the cell means.
!!
!! The mode never changes a cell mean, so the means are an honest reference, and they do not
!! depend on the bed or on the state of flotation.  Every contiguous run of 2 to 5 usable cells
!! that holds the host gives one estimate: the actual tilt less the tilt those means imply.  A
!! run of 3 cells or fewer gains one more, the tilt itself, which is the reference that assumes
!! no slope at all.  The caller reduces the set with the same minmod and the same spread that
!! the stencil readings use.
pure subroutine dg_ref_readings(tv, hv, ok, rd, nrd, lo, hi)
  real,    dimension(-4:4), intent(in)    :: tv  !< Cell tilt of the field [Z ~> m]
  real,    dimension(-4:4), intent(in)    :: hv  !< Cell mean of the field [Z ~> m]
  logical, dimension(-4:4), intent(in)    :: ok  !< True where the cell is usable
  real,    dimension(:),    intent(inout) :: rd  !< The set of estimates [Z ~> m]
  integer,                  intent(inout) :: nrd !< How many of them are set
  integer,                  intent(out)   :: lo  !< First offset of the usable run
  integer,                  intent(out)   :: hi  !< Last offset of the usable run

  real, dimension(5) :: w  ! The weights of one run [nondim]
  real :: r                ! The tilt one run's means imply [Z ~> m]
  integer :: a, b, n, k

  ! The usable run that holds the host cell.
  lo = 0 ; do k = 0, -4, -1 ; if (ok(k)) then ; lo = k ; else ; exit ; endif ; enddo
  hi = 0 ; do k = 0, 4      ; if (ok(k)) then ; hi = k ; else ; exit ; endif ; enddo
  if (hi - lo < 2) return    ! 2 cells cannot separate the mode from curvature

  do a = lo, 0 ; do b = 0, hi
    n = (b - a) + 1
    if ((n < 2) .or. (n > 5)) cycle
    w = dg_ref_w(n, -a)
    ! Summed in pairs about the middle of the run.  A quarter turn reverses the run, and a
    ! reversed running sum would not reproduce the same rounding; a symmetric pairing does.
    r = 0.0
    do k = 1, n/2
      r = r + ((w(k)*hv(a+k-1)) + (w(n+1-k)*hv(b-k+1)))
    enddo
    if (mod(n, 2) == 1) r = r + (w((n+1)/2)*hv(a + ((n-1)/2)))
    nrd = nrd + 1
    rd(nrd) = tv(0) - r
  enddo ; enddo

  ! A short run cannot see the third degree, so add the weakest reference of all as a partner
  ! for the minmod.  Without it a 3-cell run damps a smooth wave at 0.084 of its amplitude;
  ! with it, at 0.0025.
  if (hi - lo <= 2) then
    nrd = nrd + 1
    rd(nrd) = tv(0)
  endif
end subroutine dg_ref_readings

pure subroutine dg_agree_readings(gv, ok, scl, rd, nrd, cv, ncv)
  real,    dimension(-4:4), intent(in)    :: gv    !< The field per unit length [Z L-1 ~> nondim]
  logical, dimension(-4:4), intent(in)    :: ok    !< True where the cell is usable
  real,                     intent(in)    :: scl   !< Length of the centre cell along the line,
                                                   !! which returns the reading to the units of
                                                   !! the coefficient [L ~> m]
  real, dimension(:),       intent(inout) :: rd    !< The set of readings [Z ~> m]
  integer,                  intent(inout) :: nrd   !< How many of them are set
  real, dimension(:),       intent(inout) :: cv    !< The centred readings alone [Z ~> m]
  integer,                  intent(inout) :: ncv   !< How many of those are set

  real :: half_scl ! Half of that length [L ~> m]

  half_scl = 0.5*scl
  if (ok(-1) .and. ok(1)) then
    nrd = nrd + 1
    rd(nrd) = scl * (gv(0) - (0.5*(gv(-1) + gv(1))))
    ! The centred reading is kept apart because it, alone, measures the mode at this cell.
    ncv = ncv + 1
    cv(ncv) = rd(nrd)
  endif
  if (ok(-2) .and. ok(-1)) then
    nrd = nrd + 1
    rd(nrd) = half_scl * ((gv(0) - (2.0*gv(-1))) + gv(-2))
  endif
  if (ok(1) .and. ok(2)) then
    nrd = nrd + 1
    rd(nrd) = half_scl * ((gv(0) - (2.0*gv(1))) + gv(2))
  endif
end subroutine dg_agree_readings

!> The paired reading of one field. The two one-sided readings are added before they are
!! compared with the centred one. A single reading that has settled on zero then no longer
!! stops the damper, while the cancellation on smooth ice survives: for a smooth slope the
!! centred reading is -c and the one-sided ones are c-3e and c+3e, so their sum keeps the
!! sign opposite to the centred reading whatever the cubic term does.
pure subroutine dg_pair_reading(gv, ok, scl, pr, npr)
  real,    dimension(-4:4), intent(in)    :: gv  !< The field per unit length [Z L-1 ~> nondim]
  logical, dimension(-4:4), intent(in)    :: ok  !< True where the cell is usable
  real,                     intent(in)    :: scl !< Length of the centre cell, which returns the
                                                 !! reading to the units of the coefficient
                                                 !! [L ~> m]
  real, dimension(:),       intent(inout) :: pr  !< The paired reading of each field [Z ~> m]
  integer,                  intent(inout) :: npr !< How many of them are set

  real, dimension(2) :: v  ! The centred reading, then the pair [Z ~> m]
  real :: half_scl         ! Half the length of the centre cell [L ~> m]
  real :: a                ! The smaller of the two in magnitude, or zero [Z ~> m]
  real :: sp               ! Their spread, which the full set already reports [Z ~> m]
  logical :: hl, hr        ! Whether the left and the right stencil exist

  ! The pair says nothing on its own: without the centred reading there is nothing to compare
  ! it with, and without a one-sided reading there is no pair.
  if (.not.(ok(-1) .and. ok(1))) return
  hl = ok(-2) .and. ok(-1)
  hr = ok(1) .and. ok(2)
  if (.not.(hl .or. hr)) return

  half_scl = 0.5*scl
  v(1) = scl * (gv(0) - (0.5*(gv(-1) + gv(1))))
  if (hl .and. hr) then
    ! Added as one pair, so that a quarter turn, which exchanges the two, maps the sum onto
    ! itself to the last digit.
    v(2) = 0.5 * ((half_scl * ((gv(0) - (2.0*gv(-1))) + gv(-2))) + &
                  (half_scl * ((gv(0) - (2.0*gv(1))) + gv(2))))
  elseif (hl) then
    v(2) = half_scl * ((gv(0) - (2.0*gv(-1))) + gv(-2))
  else
    v(2) = half_scl * ((gv(0) - (2.0*gv(1))) + gv(2))
  endif

  call dg_minmod_spread(2, v, a, sp)
  npr = npr + 1
  pr(npr) = a
end subroutine dg_pair_reading

!> True where every usable cell of this stencil has the same grounding status. The grounded
!! fraction is set to exactly 0 and exactly 1 away from the contour, so the test needs no
!! tolerance; a cell that is part grounded counts as neither.
pure function dg_same_ground(fv, ok) result(same)
  real,    dimension(-2:2), intent(in) :: fv !< Grounded fraction [nondim]
  logical, dimension(-2:2), intent(in) :: ok !< True where the cell is usable
  logical :: same

  integer :: k, ngr, nfl

  ngr = 0 ; nfl = 0
  do k = -2, 2
    if (ok(k)) then
      if (fv(k) > 0.0) ngr = ngr + 1
      if (fv(k) < 1.0) nfl = nfl + 1
    endif
  enddo
  same = (ngr == 0) .or. (nfl == 0)
end function dg_same_ground

!> Stencil agreement for one slope direction. Read the second difference of the slope on each
!! group of three adjacent cells that contains this cell, on the thickness and again on the
!! surface form, then keep only the part every reading agrees on. A zigzag looks the same from
!! every group and gives equal readings; a real feature does not, and the readings differ. On
!! smooth ice the centred reading takes the sign opposite to the two one-sided ones, so the
!! result is exactly zero however steep the ice is.
pure subroutine dg_agree_1d(tv, bv, fv, hv, ilv, len0, ok, single_rule, edge_rule, &
                            reduce_rule, field_rule, A, S, C, nc, valid)
  real,    dimension(-4:4), intent(in)  :: tv  !< Cell tilt of the thickness [Z ~> m]
  real,    dimension(-4:4), intent(in)  :: bv  !< Cell tilt of the bed depth [Z ~> m]
  real,    dimension(-4:4), intent(in)  :: fv  !< Grounded fraction [nondim]
  real,    dimension(-4:4), intent(in)  :: hv  !< Cell mean thickness, read by the edge rule
                                               !! only [Z ~> m]
  real,    dimension(-4:4), intent(in)  :: ilv !< Reciprocal cell length along the direction,
                                               !! which the grid already carries, so the slope
                                               !! per unit length costs no division [L-1 ~> m-1]
  real,                     intent(in)  :: len0 !< Length of the centre cell [L ~> m]
  logical, dimension(-4:4), intent(in)  :: ok  !< True where the cell is usable
  logical,                  intent(in)  :: single_rule !< If true, leave a cell that has one
                                             !! stencil alone when that stencil crosses the
                                             !! grounding line
  integer,                  intent(in)  :: edge_rule !< What to do where the stencils cannot
                                             !! cancel on smooth ice; see DG1_TILT_DAMP_EDGE
  integer,                  intent(in)  :: reduce_rule !< How to reduce the readings to one
                                             !! value; see DG1_TILT_DAMP_REDUCE
  integer,                  intent(in)  :: field_rule !< Which fields to read; see
                                             !! DG1_TILT_DAMP_FIELDS
  real,                     intent(out) :: A     !< The agreed detector value [Z ~> m]
  real,                     intent(out) :: S     !< Spread of the readings [Z ~> m]
  real,                     intent(out) :: C     !< The agreed centred reading, which measures
                                                 !! the mode at this cell [Z ~> m]
  integer,                  intent(out) :: nc    !< How many centred readings were found
  logical,                  intent(out) :: valid !< True if any stencil exists

  real, dimension(16) :: rd ! Up to three stencils on two fields, or the runs of the edge rule
                            ! [Z ~> m]
  real, dimension(2) :: cv  ! The centred reading of each field [Z ~> m]
  real, dimension(2) :: pr  ! The paired reading of each field [Z ~> m]
  real, dimension(-4:4) :: gv ! One field per unit length [Z L-1 ~> nondim]
  real :: Sc                ! Spread of the centred readings, not used [Z ~> m]
  real :: Sp                ! Spread of the paired readings, not used [Z ~> m]
  integer :: k, nrd, nst, lo, hi, npr
  logical :: st_c, st_l, st_r ! Whether the centred, left and right stencils exist
  logical :: paired         ! Whether to add the one-sided readings before comparing them
  logical :: use_t, use_s   ! Which of the two fields this cell reads

  A = 0.0 ; S = 0.0 ; C = 0.0 ; nc = 0 ; valid = .false.
  st_c = ok(-1) .and. ok(1)
  st_l = ok(-2) .and. ok(-1)
  st_r = ok(1) .and. ok(2)
  nst = 0
  if (st_c) nst = nst + 1
  if (st_l) nst = nst + 1
  if (st_r) nst = nst + 1

  ! The cancellation on smooth ice works through the opposite signs of the centred reading and
  ! the one-sided ones, so it needs the centred stencil AND at least one of the others.  Where
  ! it is not available the cell means take over.
  if ((edge_rule > 0) .and. .not.(st_c .and. (st_l .or. st_r))) then
    nrd = 0
    call dg_ref_readings(tv, hv, ok, rd, nrd, lo, hi)
    if (nrd == 0) return
    ! The means themselves bend where the run crosses a flotation break, so the reference is
    ! not to be trusted there.  The grounded fraction is set to exactly 0 and exactly 1, so
    ! the test needs no tolerance.
    ! DAMP_EDGE_REF_GL reads the means across the break as well.  The means themselves do
    ! bend there, but leaving the cell alone is not free either: across a channel the run is
    ! the whole column, so one partly grounded cell silences both walls.
    if (edge_rule == DAMP_EDGE_REF) then
      if (maxval(fv(lo:hi)) > minval(fv(lo:hi))) return
    endif
    valid = .true.
    call dg_minmod_spread(nrd, rd, A, S)
    return
  endif

  if (nst == 0) return
  valid = .true.

  ! One stencil is no second opinion: the readings of the two fields are all there is, and on
  ! floating ice they are the same reading twice.  Where that stencil also crosses the grounding
  ! line it reads the change of slope there, so leave the cell alone.  The test is exact because
  ! the grounded fraction is set to exactly 0 and exactly 1 away from the contour.
  if (single_rule .and. (nst == 1)) then
    if (st_l) then
      if (.not.(all(fv(-2:0) <= 0.0) .or. all(fv(-2:0) >= 1.0))) return
    elseif (st_r) then
      if (.not.(all(fv(0:2) <= 0.0) .or. all(fv(0:2) >= 1.0))) return
    else
      if (.not.(all(fv(-1:1) <= 0.0) .or. all(fv(-1:1) >= 1.0))) return
    endif
  endif

  ! Add the one-sided readings only where the cells this stencil reads share a grounding
  ! status.  Where they do not, the surface has a real kink, and one kinked cell reads exactly
  ! like the mode; there the comparison of every reading with every other one is what leaves
  ! the grounding line in place.
  paired = .false.
  if (reduce_rule == DAMP_RED_PAIRED) paired = dg_same_ground(fv(-2:2), ok(-2:2))

  nrd = 0 ; nc = 0 ; npr = 0
  ! The thickness, then the surface form: the thickness tilt less the part the bed explains,
  ! which is the change of surface elevation across a grounded cell and the tilt itself afloat.
  ! Where the ice floats the bed term is zero, so the two fields are one field and the
  ! selection changes nothing.  It bites on grounded ice only.
  use_t = (field_rule /= DAMP_FLD_SURF)
  use_s = (field_rule /= DAMP_FLD_THCK)
  if (field_rule == DAMP_FLD_GL_TH) then
    ! A cell the contour crosses has a grounded fraction strictly between 0 and 1, so its
    ! surface form removes a FRACTION of the bed.  That is no surface, and holding the
    ! thickness against it only silences the cell.  Drop it and let the thickness decide.
    use_t = .true.
    use_s = .not.((fv(0) > 0.0) .and. (fv(0) < 1.0))
  endif
  if (use_t) then
    do k = -2, 2
      gv(k) = tv(k) * ilv(k)
    enddo
    call dg_agree_readings(gv, ok, len0, rd, nrd, cv, nc)
    if (paired) call dg_pair_reading(gv, ok, len0, pr, npr)
  endif
  if (use_s) then
    do k = -2, 2
      gv(k) = (tv(k) - (fv(k)*bv(k))) * ilv(k)
    enddo
    call dg_agree_readings(gv, ok, len0, rd, nrd, cv, nc)
    if (paired) call dg_pair_reading(gv, ok, len0, pr, npr)
  endif

  ! The spread always comes from the full set, because it reports how far the readings stand
  ! apart, which the pair hides by construction.
  call dg_minmod_spread(nrd, rd, A, S)
  ! The two fields must still agree: a thickness that follows a rough bed is protected here
  ! exactly as it is in the full set.
  if (paired .and. (npr > 0)) call dg_minmod_spread(npr, pr, A, Sp)
  ! The two fields must still agree on the centred reading, so that a thickness which follows
  ! a rough bed is protected exactly as it is in the full set.
  if (nc > 0) call dg_minmod_spread(nc, cv, C, Sc)
end subroutine dg_agree_1d

!> Stencil agreement for the twist. The checkerboard alternates along both axes, so it is read
!! along both, and it is constant along either diagonal, so a diagonal stencil would read zero
!! and silence the damper. The field per unit length is the twist per unit area; along a line of
!! constant j only the widths vary, so each line reduces to the same one-dimensional form.
pure subroutine dg_agree_2d(wv, bwv, fv, ilx, ily, dx0, dy0, okw, edge_rule, reduce_rule, &
                            field_rule, A, S, C, nc, valid)
  real,    dimension(-4:4,-4:4), intent(in) :: wv  !< Cell twist of the thickness [Z ~> m]
  real,    dimension(-4:4,-4:4), intent(in) :: bwv !< Cell twist of the bed depth [Z ~> m]
  real,    dimension(-4:4,-4:4), intent(in) :: fv  !< Grounded fraction [nondim]
  real,    dimension(-4:4), intent(in) :: ilx !< Reciprocal cell width along the row [L-1 ~> m-1]
  real,    dimension(-4:4), intent(in) :: ily !< Reciprocal cell height along the column
                                              !! [L-1 ~> m-1]
  real,                     intent(in) :: dx0 !< Width of the centre cell [L ~> m]
  real,                     intent(in) :: dy0 !< Height of the centre cell [L ~> m]
  logical, dimension(-4:4,-4:4), intent(in) :: okw !< True where the cell is usable
  integer,                  intent(in)  :: edge_rule !< What to do where neither axis can
                                             !! cancel on smooth ice; see DG1_TILT_DAMP_EDGE
  integer,                  intent(in)  :: reduce_rule !< How to reduce the readings to one
                                             !! value; see DG1_TILT_DAMP_REDUCE
  integer,                  intent(in)  :: field_rule !< Which fields to read; see
                                             !! DG1_TILT_DAMP_FIELDS
  real,                     intent(out) :: A     !< The agreed detector value [Z ~> m]
  real,                     intent(out) :: S     !< Spread of the readings [Z ~> m]
  real,                     intent(out) :: C     !< The agreed centred reading, which measures
                                                 !! the mode at this cell [Z ~> m]
  integer,                  intent(out) :: nc    !< How many centred readings were found
  logical,                  intent(out) :: valid !< True if any stencil exists

  real, dimension(12) :: rd ! Up to six stencils on two fields [Z ~> m]
  real, dimension(4) :: cv  ! The centred reading of each field on each axis [Z ~> m]
  real, dimension(4) :: pr  ! The paired reading of each field on each axis [Z ~> m]
  real, dimension(-4:4) :: gv ! One field per unit length along the line [Z L-1 ~> nondim]
  real, dimension(-2:2) :: fl ! Grounded fraction along the line [nondim]
  logical, dimension(-4:4) :: okl ! Usability along the line
  real :: Sc                ! Spread of the centred readings, not used [Z ~> m]
  real :: Sp                ! Spread of the paired readings, not used [Z ~> m]
  integer :: k, nrd, npr
  logical :: st_x, st_y     ! Whether each axis can cancel on smooth ice
  logical :: pr_x, pr_y     ! Whether each axis may add its one-sided readings
  logical :: use_t, use_s   ! Which of the two fields this cell reads

  A = 0.0 ; S = 0.0 ; C = 0.0 ; nc = 0 ; valid = .false.
  nrd = 0 ; npr = 0
  use_t = (field_rule /= DAMP_FLD_SURF)
  use_s = (field_rule /= DAMP_FLD_THCK)
  if (field_rule == DAMP_FLD_GL_TH) then
    use_t = .true.
    use_s = .not.((fv(0,0) > 0.0) .and. (fv(0,0) < 1.0))
  endif

  ! The twist puts both axes into one set, so one axis that can cancel protects the whole set.
  ! Only a corner loses both.  There the cell means would have to supply a twist, which needs a
  ! two-dimensional fit that has not been tested, so stand down instead.
  st_x = (okw(-1,0) .and. okw(1,0)) .and. &
         ((okw(-2,0) .and. okw(-1,0)) .or. (okw(1,0) .and. okw(2,0)))
  st_y = (okw(0,-1) .and. okw(0,1)) .and. &
         ((okw(0,-2) .and. okw(0,-1)) .or. (okw(0,1) .and. okw(0,2)))
  if ((edge_rule > 0) .and. .not.(st_x .or. st_y)) return

  ! Along x, at the row of this cell.  The twist per unit area is w/(dx*dy), and dy does not
  ! change along a row, so the height divides out of the second difference and the scale alike.
  do k = -2, 2
    okl(k) = okw(k,0)
    gv(k) = wv(k,0) * ilx(k)
  enddo
  ! Each axis decides for itself: a row that crosses the grounding line must not add its
  ! one-sided readings even where the column may.
  pr_x = .false.
  if (reduce_rule == DAMP_RED_PAIRED) then
    do k = -2, 2
      fl(k) = fv(k,0)
    enddo
    pr_x = dg_same_ground(fl, okl(-2:2))
  endif
  if (use_t) then
    call dg_agree_readings(gv, okl, dx0, rd, nrd, cv, nc)
    if (pr_x) call dg_pair_reading(gv, okl, dx0, pr, npr)
  endif
  if (use_s) then
    do k = -2, 2
      gv(k) = (wv(k,0) - (fv(k,0)*bwv(k,0))) * ilx(k)
    enddo
    call dg_agree_readings(gv, okl, dx0, rd, nrd, cv, nc)
    if (pr_x) call dg_pair_reading(gv, okl, dx0, pr, npr)
  endif

  ! Along y, at the column of this cell.
  do k = -2, 2
    okl(k) = okw(0,k)
    gv(k) = wv(0,k) * ily(k)
  enddo
  pr_y = .false.
  if (reduce_rule == DAMP_RED_PAIRED) then
    do k = -2, 2
      fl(k) = fv(0,k)
    enddo
    pr_y = dg_same_ground(fl, okl(-2:2))
  endif
  if (use_t) then
    call dg_agree_readings(gv, okl, dy0, rd, nrd, cv, nc)
    if (pr_y) call dg_pair_reading(gv, okl, dy0, pr, npr)
  endif
  if (use_s) then
    do k = -2, 2
      gv(k) = (wv(0,k) - (fv(0,k)*bwv(0,k))) * ily(k)
    enddo
    call dg_agree_readings(gv, okl, dy0, rd, nrd, cv, nc)
    if (pr_y) call dg_pair_reading(gv, okl, dy0, pr, npr)
  endif

  valid = (nrd > 0)
  ! The spread always comes from the full set; the pair hides it by construction.
  call dg_minmod_spread(nrd, rd, A, S)
  if ((npr > 0) .and. (reduce_rule == DAMP_RED_PAIRED)) call dg_minmod_spread(npr, pr, A, Sp)
  ! Both axes and both fields must agree on the centred reading, exactly as they must on the
  ! full set.  A minmod over a set does not depend on the order a quarter turn presents it in.
  if (nc > 0) call dg_minmod_spread(nc, cv, C, Sc)
end subroutine dg_agree_2d

!> Nodal rates damping the grid-scale in-cell tilt and twist. A tilt alternating between cells
!! adds nothing to any face jump, so the upwind flux and the artificial viscosity cannot see it,
!! but it gives a spurious surface slope. The detector is the tilt Laplacian
!! A = t_j - (t_j-1 + t_j+1)/2, zero on a uniform tilt and O(dx^3) when smooth. Only the part
!! not explained by the bed (grounded share only) or by the tilt of the cell means is removed,
!! gated by ds/dh relative to the thickness. Grounding-line cells are protected per
!! DG1_TILT_DAMP_GL_REACH unless the kink reference explains them. The twist
!! (SW - SE - NW + NE) is treated the same way with a 2D detector. The corrections are
!! equal and opposite between nodes, so the cell mean is unchanged.
subroutine dg_nodal_mode_damp_rate(CS, G, hmask, h_nodal_in, dt, T_node)
  type(ice_shelf_dyn_CS), intent(in)  :: CS   !< Ice shelf dynamics control structure
  type(ocean_grid_type),  intent(in)  :: G    !< Ocean grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in) :: hmask !< Cell mask
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(in)  :: h_nodal_in !< Corner thicknesses [Z ~> m]
  real,                   intent(in)  :: dt   !< Time step [T ~> s]
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(out) :: T_node !< Nodal tilt-damping rate [Z T-1 ~> m s-1]

  real, dimension(SZDI_(G),SZDJ_(G)) :: t_xi   ! Per-cell xi-tilt of thickness [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: t_eta  ! Per-cell eta-tilt of thickness [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: b_xi   ! Per-cell xi-tilt of the bed [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: b_eta  ! Per-cell eta-tilt of the bed [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: hbar_c ! Mean of the four corner thicknesses [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: sbar_c ! Surface elevation from the cell mean thickness,
                                  ! the reference of the surface-slope gate [Z ~> m]
  logical, dimension(SZDI_(G),SZDJ_(G)) :: ice_ok !< True on a fully ice-covered cell
  real, dimension(SZDI_(G),SZDJ_(G)) :: w_c    ! Per-cell xy-twist of thickness [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: b_w    ! Per-cell xy-twist of the bed [Z ~> m]
  real :: A_h      ! Tilt Laplacian of the thickness [Z ~> m]
  real :: A_bed    ! Tilt Laplacian of the bed, zero on a floating cell [Z ~> m]
  real :: A_ref    ! Tilt Laplacian of the mean-supported tilt [Z ~> m]
  real :: A_w      ! 4-neighbour Laplacian of the twist [Z ~> m]
  real :: A_w_bed  ! As A_w for the bed [Z ~> m]
  real :: A_w_ref  ! As A_w for the mean-supported twist [Z ~> m]
  real :: href_d   ! Gate-normalizing thickness for the direction in hand [Z ~> m]
  real :: href_w   ! Gate-normalizing thickness for the twist stencil [Z ~> m]
  real, dimension(-4:4) :: tv, bv, hv ! Stencil gathers of tilt, bed tilt and mean
  real, dimension(-4:4) :: gv         ! Stencil gather of the grounded fraction [nondim]
  logical, dimension(-4:4) :: okv     ! Stencil gather of usability
  logical :: det_ok ! True if the detector found an admissible stencil
  real :: kink_c ! How far the reference carried the flotation break itself [nondim]
  integer :: dmode  ! Stencil the detector chose: 1 centred, 2 forward, 3 backward
  real :: gwt       ! Grounding-line weight over the reach of that stencil [nondim]
  real :: A_dmp     ! The part of the detector the correction actually removes [Z ~> m]
  real :: kqx, kqy, kq0 ! The flotation break near the cell as one linear function [Z ~> m]
  real :: ksig      ! Span of the grounded fraction over the twist's window [nondim]
  real :: kconf     ! How far a single-line fit describes this neighbourhood [nondim]
  logical :: kw_ok  ! True if that break could be reconstructed at all
  real :: knat_xi, knat_eta ! Removal rate the transport already delivers on each axis [T-1 ~> s-1]
  real :: kwant     ! Rate the gate asks for, before the transport's share [T-1 ~> s-1]
  real :: itau_xi, itau_eta, itau_w ! Full rate per mode, from DG1_TILT_DAMP_U_CUT [T-1 ~> s-1]
  real :: d_ahat    ! Sum of squared detector outputs, for diagnostics [Z2 ~> m2]
  real :: d_spread  ! Largest reading spread over the three modes, for diagnostics [Z ~> m]
  real :: d_gate    ! Largest gate over the three modes [nondim]
  real :: d_want, d_got ! Rates summed over modes, before and after the transport credit [T-1 ~> s-1]
  logical :: diag_on ! True if any of the three activity diagnostics is registered
  integer :: k, kk, kj, jj ! Stencil offsets, and the clamped array indices
  ! Twist gathers, over the 7x7 block the biased cross difference can reach.
  real, dimension(-4:4,-4:4) :: ww, bw, hw, fw ! Twist, bed twist, cell mean, grounded frac
  real, dimension(-4:4) :: gwx, gwy        ! Grounded fraction along each axis [nondim]
  logical, dimension(-4:4,-4:4) :: okw     ! Usability
  real, dimension(-4:4) :: wrx, wry ! Mean-supported twist along each axis [Z ~> m]
  integer :: bx, by, m              ! Bias along xi, along eta, and the trial index
  logical :: tw_ok                  ! True once an admissible bias pair is found
  integer :: ntw, itw, tier         ! Admissible pairs in the tier in use, and loop indices
  integer, dimension(4) :: bx_use, by_use ! Those pairs
  real :: reach_w                   ! Grounding-line reach weight, worst over those pairs
  real :: inv_ntw                   ! Reciprocal of the number of pairs [nondim]
  ! Offsets a bias reaches: S for the second difference, F for the first.
  integer, dimension(3), parameter :: Slo = (/ -2, 0, -4 /), Shi = (/ 2, 4, 0 /)
  integer, dimension(3), parameter :: Flo = (/ -1, 0, -2 /), Fhi = (/  1, 2, 0 /)
  ! The two NON-centre cells the second difference uses, for the stencil mean.
  integer, dimension(3), parameter :: Elo = (/ -1, 1, -1 /), Ehi = (/  1, 2, -2 /)
  ! Trial order: centred in both axes first, then centred in one, then neither.
  integer, dimension(9), parameter :: bx_try = (/ 1,1,1, 2,3, 2,2,3,3 /)
  integer, dimension(9), parameter :: by_try = (/ 1,2,3, 1,1, 2,3,2,3 /)
  ! Tiers of that list: centred in both, biased in one, biased in both.
  integer, dimension(3), parameter :: tw_lo = (/ 1, 2, 6 /), tw_hi = (/ 1, 5, 9 /)
  real :: floor_A  ! Larger of the bed- and mean-supported floors [Z ~> m]
  real :: excess   ! One-sided excess over that floor [Z ~> m]
  real :: href     ! Cell-mean thickness used to normalize the gate [Z ~> m]
  real :: dsdh     ! ds/dh, 1 grounded and 1 - rho_i/rho_w floating [nondim]
  real :: dsurf    ! Real rise of the surface across the cell, the surface-slope gate's
                   ! reference [Z ~> m]
  real :: A_spread ! Largest minus smallest detector reading; zero unless the detector
                   ! reports it [Z ~> m]
  real :: A_cen    ! The agreed centred reading, which measures the mode at the cell [Z ~> m]
  integer :: n_cen ! How many centred readings the detector found
  real :: dx_cell, dy_cell ! Cell width and height [L ~> m]
  real, dimension(-4:4) :: ilv ! Reciprocal cell length along the direction examined [L-1 ~> m-1]
  real, dimension(-4:4) :: iwx ! Reciprocal cell width along the twist's row [L-1 ~> m-1]
  real, dimension(-4:4) :: iwy ! Reciprocal cell height along the twist's column [L-1 ~> m-1]
  real, dimension(SZDI_(G),SZDJ_(G)) :: r_xi, r_eta, r_w ! Removal per mode, before the
                                  ! optional high pass [Z T-1 ~> m s-1]
  real :: rf_xi, rf_eta, rf_w ! The same after it [Z T-1 ~> m s-1]
  integer :: isc_w, iec_w, jsc_w, jec_w ! Loop bounds, one cell wider with the high pass
  logical :: skip_xi, skip_eta, skip_w ! True where the transport already delivers the full rate
  integer :: kglo, kghi ! Gather range: the agreement detector reads two cells, not four
  logical :: agree_det  ! True when the stencil-agreement detector is selected
  logical :: hybrid     ! True when every cell chooses its own detector
  logical :: agr_xi, agr_eta, agr_w ! True where this component uses stencil agreement
  integer :: gate_now   ! The gate that belongs to the detector this component used
  integer :: ngr, nfl   ! Cells of the twist window that are not afloat, and not grounded
  integer :: ki, kj_g   ! Offsets across the twist window
  logical :: edge_on    ! True when the edge rule reads the cell means
  logical :: slope_gate ! True when the gate reads the real surface slope
  real :: gam      ! Gate fraction, 0 to 1 [nondim]
  real :: kap      ! Delivered rate for this cell and direction [T-1 ~> s-1]
  real :: rate_cap ! Largest rate the explicit step may carry [T-1 ~> s-1]
  real :: c_ucut   ! Rate the transport cut-off speed asks for, per unit of reciprocal
                   ! cell length [L T-1 ~> m s-1]
  real :: rhoi_rhow ! Ice/ocean density ratio for the flotation test [nondim]
  real :: one_m_r  ! 1 - rho_i/rho_w, the floating-ice ds/dh [nondim]
  real :: Tbar     ! Cell mean of the nodal increment, removed so the damper moves no mass
                   ! [Z T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G)) :: f_gnd ! CS%ground_frac clipped to [0,1] [nondim]
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed

  ! The caller reads the compute domain only, and skips the cells this routine skips, so there
  ! is no need to clear the halo.
  T_node(G%isc:G%iec,G%jsc:G%jec,:,:) = 0.0
  diag_on = (CS%id_dg_damp_tend > 0) .or. (CS%id_dg_damp_gate > 0) .or. &
            (CS%id_dg_damp_want > 0) .or. (CS%id_dg_damp_got > 0)
  if (diag_on) then
    CS%dg_damp_tend(:,:) = 0.0 ; CS%dg_damp_gate(:,:) = 0.0
    CS%dg_damp_want(:,:) = 0.0 ; CS%dg_damp_got(:,:) = 0.0
  endif
  if (.not.(CS%dg_tilt_damp .or. CS%dg_twist_damp)) return

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  ! A_ref reads cell means two cells away.
  if ((G%isc - G%isd < 2) .or. (G%jsc - G%jsd < 2)) call MOM_error(FATAL, &
    "dg_nodal_mode_damp_rate: DG1_TILT_DAMP needs a halo of at least 2 cells; "//&
    "increase NIHALO/NJHALO.")

  if (.not.associated(CS%bed_node)) call MOM_error(FATAL, &
    "dg_nodal_mode_damp_rate: DG1_TILT_DAMP requires a nodal bed (CS%bed_node), "//&
    "which sets both the flotation branch and the bed-supported tilt floor.")
  rhoi_rhow = CS%density_ice / CS%density_ocean_avg
  one_m_r = 1.0 - rhoi_rhow

  ! Hold the delivered rate to half the SSP-RK2 single-cell budget so this term
  ! cannot on its own drive the step past the stable range; the remainder is
  ! left to advection and the artificial viscosity.
  rate_cap = 0.5 / max(dt, tiny(dt))

  ! The full rate per unit of reciprocal cell length, which no cell changes.
  c_ucut = CS%dg_damp_advective_c * CS%dg_damp_u_cut

  agree_det  = (CS%dg_damp_detector == DAMP_DET_AGREE)
  ! The hybrid gathers as the tilt Laplacian does, because either detector may run at any
  ! cell, and it decides between them cell by cell inside the loop.
  hybrid     = (CS%dg_damp_detector == DAMP_DET_HYBRID)
  ! The edge rule fits the cell means over runs of up to five cells, so it reaches as far as
  ! today's detector does and it needs the means that the agreement path otherwise skips.
  edge_on    = (agree_det .or. hybrid) .and. (CS%dg_damp_edge_rule > 0)
  slope_gate = (CS%dg_damp_gate_form == DAMP_GATE_SLOPE)

  ! Cell-local fields on the full data domain.  The loop below writes every element of
  ! hbar_c, ice_ok, t_xi, t_eta and f_gnd, so clearing them first would only be a dead
  ! store; w_c and sbar_c are written only where their feature asks for them, and nothing
  ! reads them otherwise.
  do j = jsd, jed ; do i = isd, ied
    hbar_c(i,j) = 0.25*((h_nodal_in(i,j,1,1) + h_nodal_in(i,j,2,2)) + &
                        (h_nodal_in(i,j,2,1) + h_nodal_in(i,j,1,2)))
    ice_ok(i,j) = (hmask(i,j) == 1.0)
    if (CS%dg_tilt_damp) then
      t_xi(i,j)  = 0.5*((h_nodal_in(i,j,2,1) - h_nodal_in(i,j,1,1)) + &
                        (h_nodal_in(i,j,2,2) - h_nodal_in(i,j,1,2)))
      t_eta(i,j) = 0.5*((h_nodal_in(i,j,1,2) - h_nodal_in(i,j,1,1)) + &
                        (h_nodal_in(i,j,2,2) - h_nodal_in(i,j,2,1)))
    endif
    if (CS%dg_twist_damp) &
      w_c(i,j) = (h_nodal_in(i,j,2,2) - h_nodal_in(i,j,1,2)) - &
                 (h_nodal_in(i,j,2,1) - h_nodal_in(i,j,1,1))
    ! The sub-element grounded fraction
    ! TODO: use f_ground_cell instead for CISM-style grounding with DG for advection?
    f_gnd(i,j) = min(max(CS%ground_frac(i,j), 0.0), 1.0)
    ! Surface elevation the cell mean implies: freeboard afloat, thickness above the bed where
    ! the flotation deficit rho_i/rho_w*h - depth is positive. Continuous across the contour.
    ! Only the surface-slope gate reads it.
    if (slope_gate) &
      sbar_c(i,j) = (one_m_r*hbar_c(i,j)) + max((rhoi_rhow*hbar_c(i,j)) - CS%bed_elev(i,j), 0.0)
  enddo ; enddo

  ! Bed tilts read bed_node, so they stop one ring short; no chosen stencil reads that ring,
  ! but the gathers clamp into it, so clear the ring alone rather than the whole array.
  if (CS%dg_tilt_damp) then
    b_xi(isd,:) = 0.0 ; b_xi(ied,:) = 0.0 ; b_xi(:,jsd) = 0.0 ; b_xi(:,jed) = 0.0
    b_eta(isd,:) = 0.0 ; b_eta(ied,:) = 0.0 ; b_eta(:,jsd) = 0.0 ; b_eta(:,jed) = 0.0
  endif
  if (CS%dg_twist_damp) then
    b_w(isd,:) = 0.0 ; b_w(ied,:) = 0.0 ; b_w(:,jsd) = 0.0 ; b_w(:,jed) = 0.0
  endif
  do j = jsd+1, jed-1 ; do i = isd+1, ied-1
    if (CS%dg_tilt_damp) then
      b_xi(i,j)  = 0.5*((CS%bed_node(I,J-1) - CS%bed_node(I-1,J-1)) + &
                        (CS%bed_node(I,J)   - CS%bed_node(I-1,J)))
      b_eta(i,j) = 0.5*((CS%bed_node(I-1,J) - CS%bed_node(I-1,J-1)) + &
                        (CS%bed_node(I,J)   - CS%bed_node(I,J-1)))
    endif
    if (CS%dg_twist_damp) &
      b_w(i,j) = (CS%bed_node(I,J) - CS%bed_node(I-1,J)) - &
                 (CS%bed_node(I,J-1) - CS%bed_node(I-1,J-1))
  enddo ; enddo

  ! The high pass reads the removal of the neighbours, so the removal is built one cell beyond
  ! the compute domain, and the stencils behind it reach two cells further still.  The bed tilts
  ! stop one ring inside the data domain, which is what sets the halo this needs.
  isc_w = isc ; iec_w = iec ; jsc_w = jsc ; jec_w = jec
  if (CS%dg_damp_filter) then
    if ((G%isc - G%isd < 4) .or. (G%jsc - G%jsd < 4)) call MOM_error(FATAL, &
      "dg_nodal_mode_damp_rate: DG1_TILT_DAMP_FILTER needs a halo of at least 4 cells; "//&
      "increase NIHALO/NJHALO or switch the filter off.")
    isc_w = isc-1 ; iec_w = iec+1 ; jsc_w = jsc-1 ; jec_w = jec+1
  endif
  ! Only the window the second pass reads has to be valid.
  r_xi(isc_w:iec_w,jsc_w:jec_w) = 0.0
  r_eta(isc_w:iec_w,jsc_w:jec_w) = 0.0
  r_w(isc_w:iec_w,jsc_w:jec_w) = 0.0
  kglo = -4 ; kghi = 4
  if (agree_det .and. (.not.edge_on)) then ; kglo = -2 ; kghi = 2 ; endif

  do j = jsc_w, jec_w ; do i = isc_w, iec_w
    if (hmask(i,j) /= 1.0) cycle
    href = hbar_c(i,j)
    if (href <= 0.0) cycle
    if (diag_on) then
      d_gate = 0.0 ; d_want = 0.0 ; d_got = 0.0
      d_ahat = 0.0 ; d_spread = 0.0
    endif
    ! The grid carries both the length and its reciprocal, so neither costs a division.
    dx_cell = G%dxT(i,j) ; dy_cell = G%dyT(i,j)
    ! Only a detector that reports the spread of its readings can open the agreement gate on
    ! anything but a perfect zigzag; today's detector leaves it zero.
    A_spread = 0.0

    ! The zigzag in the thickness raises the surface by less than itself where the ice
    ! floats.  Both the thickness gate and the slope gate compare surface rises; the
    ! agreement gate compares readings with each other and never reads this.
    dsdh = 0.0
    if (CS%dg_damp_gate_form /= DAMP_GATE_AGREEMENT) &
      dsdh = one_m_r + f_gnd(i,j)*(1.0 - one_m_r)

    ! Full rate c*U_CUT/dx along each axis; the twist uses sqrt(dx*dy). knat is the rate the
    ! transport itself delivers, from the cell-centred |u| or |v|.
    itau_xi  = c_ucut * G%IdxT(i,j)
    itau_eta = c_ucut * G%IdyT(i,j)
    itau_w   = 0.0
    if (CS%dg_twist_damp) itau_w = c_ucut * sqrt(G%IdxT(i,j) * G%IdyT(i,j))

    knat_xi = 0.0 ; knat_eta = 0.0
    if (CS%dg_damp_advective) then
      knat_xi = CS%dg_damp_advective_c * G%IdxT(i,j) * 0.25 * &
                ((abs(CS%u_shelf(I-1,J-1)) + abs(CS%u_shelf(I,J))) + &
                 (abs(CS%u_shelf(I,J-1))   + abs(CS%u_shelf(I-1,J))))
      knat_eta = CS%dg_damp_advective_c * G%IdyT(i,j) * 0.25 * &
                ((abs(CS%v_shelf(I-1,J-1)) + abs(CS%v_shelf(I,J))) + &
                 (abs(CS%v_shelf(I,J-1))   + abs(CS%v_shelf(I-1,J))))
    endif

    ! The delivered rate is min(max(gam*itau - knat, 0), cap) with gam <= 1, so where the
    ! transport already removes the mode faster than the damper's full rate there is nothing to
    ! add and the detector need not run at all.  Guarded on the diagnostics, which report the
    ! rate the gate asked for and would otherwise lose those cells.
    skip_xi  = (.not.diag_on) .and. (itau_xi  <= knat_xi)
    skip_eta = (.not.diag_on) .and. (itau_eta <= knat_eta)
    skip_w   = (.not.diag_on) .and. (itau_w   <= min(knat_xi, knat_eta))

    ! --- xi direction, gated only on row j ---
    if (CS%dg_tilt_damp .and. (.not.skip_xi)) then
      okv(-4:kglo-1) = .false. ; okv(kghi+1:4) = .false.
      do k = kglo, kghi
        kk = min(max(i+k, isd), ied)
        okv(k) = ice_ok(kk,j) .and. (i+k >= isd) .and. (i+k <= ied)
        tv(k) = t_xi(kk,j) ; bv(k) = b_xi(kk,j)
        gv(k) = f_gnd(kk,j) ; ilv(k) = G%IdxT(kk,j)
        ! Agreement needs no stencil thickness unless the edge rule reads the means.
        if (agree_det .and. (.not.edge_on)) cycle
        hv(k) = hbar_c(kk,j)
      enddo
      ! The tilt Laplacian is the stronger detector, but its reference averages the cell
      ! means across whatever the stencil covers, so a flotation break inside the stencil
      ! bends it.  Hand those cells to stencil agreement, which holds the readings against
      ! each other and needs no reference.  The grounded fraction is set to exactly 0 and
      ! exactly 1 away from the contour, so the test needs no tolerance.
      agr_xi = agree_det
      if (hybrid) then
        ! Try the tilt Laplacian first and keep the cell unless it FAILS.  Failure is
        ! either that it finds no stencil at all - it reads three cells, but the reference
        ! it is held against reads five - or that the reach exemption refuses to trust the
        ! reference it did build, which is what happens where the window straddles the
        ! flotation contour.  Both of those leave the cell undamped.  Stencil agreement
        ! needs no reference and no run of five, so it takes exactly those cells, and the
        ! stronger detector keeps every other one.
        call dg_tilt_detector_1d(tv, bv, hv, okv, gv, CS%dg_damp_kink_ref, &
                                 CS%dg_damp_kink_fit_tol, &
                                 A_h, A_bed, A_ref, href_d, det_ok, dmode, kink_c)
        agr_xi = .true.
        if (det_ok) then
          gwt = kink_c + (1.0 - kink_c) * &
                dg_gl_reach_wt(gv, max(dmode,1), CS%dg_damp_gl_reach)
          agr_xi = (gwt <= 0.0)
        endif
      endif
      if (agr_xi) then
        ! Agreement needs no floor and no kink reference, and nothing here protects the
        ! grounding line: the readings disagree there on their own.  No single stencil is
        ! chosen, so the thickness gate normalizes on the cell mean.
        call dg_agree_1d(tv, bv, gv, hv, ilv, dx_cell, okv, CS%dg_damp_single_rule, &
                         CS%dg_damp_edge_rule, CS%dg_damp_reduce, CS%dg_damp_fields, &
                         A_dmp, A_spread, A_cen, n_cen, det_ok)
        ! The minmod keeps the veto; the centred reading sets the amount.  Where no centred
        ! reading exists the cell has no second opinion at all, so leave it as it was.
        if (CS%dg_damp_centred_amp .and. (A_dmp /= 0.0) .and. (n_cen > 0)) A_dmp = A_cen
        excess = abs(A_dmp) ; gwt = 1.0 ; href_d = href
      else
        ! Under the hybrid the detector ran above and gwt is already in hand; its answer
        ! is what chose this branch.
        if (.not.hybrid) then
          call dg_tilt_detector_1d(tv, bv, hv, okv, gv, CS%dg_damp_kink_ref, &
                                   CS%dg_damp_kink_fit_tol, &
                                   A_h, A_bed, A_ref, href_d, det_ok, dmode, kink_c)
          gwt = kink_c + (1.0 - kink_c) * &
                dg_gl_reach_wt(gv, max(dmode,1), CS%dg_damp_gl_reach)
        endif
        ! One-sided, so ice carrying LESS structure than the bed forces is left
        ! alone rather than driven further from it.  max(), not a sum: over a
        ! rough bed the means already contain the bed's own structure.
        A_bed = f_gnd(i,j)*A_bed
        if (CS%dg_damp_excess_only) then
          A_dmp = dg_unexplained(A_h, A_ref, A_bed)
          excess = abs(A_dmp)
        else
          floor_A = max(abs(A_bed), abs(A_ref))
          excess = abs(A_h) - floor_A
          A_dmp = A_h
        endif
      endif
      if (det_ok .and. (excess > 0.0)) then
        dsurf = 0.0
        if (slope_gate) dsurf = 0.5*(sbar_c(i+1,j) - sbar_c(i-1,j))
        ! Each branch is gated by the gate that can read it: only the agreement detector
        ! reports a spread, and only the tilt Laplacian is normalized by a thickness.
        gate_now = CS%dg_damp_gate_form
        if (hybrid .and. agr_xi) gate_now = DAMP_GATE_AGREEMENT
        gam = gwt * dg_damp_gate_frac(gate_now, dsdh, excess, A_spread, href_d, &
                                      dsurf, CS%dg_tilt_damp_r_hi, CS%dg_damp_rho_s, &
                                      CS%dg_damp_slope_floor*dx_cell, CS%dg_damp_rho_g)
        kwant = gam * itau_xi
        kap = min(max(kwant - knat_xi, 0.0), rate_cap)
        r_xi(i,j) = kap*A_dmp
        if (diag_on) then
          d_ahat = d_ahat + (A_dmp**2) ; d_spread = max(d_spread, A_spread)
          d_gate = max(d_gate, gam) ; d_want = d_want + kwant
          d_got = d_got + kap
        endif
      endif
    endif

    ! --- eta direction ---
    if (CS%dg_tilt_damp .and. (.not.skip_eta)) then
      okv(-4:kglo-1) = .false. ; okv(kghi+1:4) = .false.
      do k = kglo, kghi
        kk = min(max(j+k, jsd), jed)
        okv(k) = ice_ok(i,kk) .and. (j+k >= jsd) .and. (j+k <= jed)
        tv(k) = t_eta(i,kk) ; bv(k) = b_eta(i,kk)
        gv(k) = f_gnd(i,kk) ; ilv(k) = G%IdyT(i,kk)
        if (agree_det .and. (.not.edge_on)) cycle
        hv(k) = hbar_c(i,kk)
      enddo
      ! The tilt Laplacian is the stronger detector, but its reference averages the cell
      ! means across whatever the stencil covers, so a flotation break inside the stencil
      ! bends it.  Hand those cells to stencil agreement, which holds the readings against
      ! each other and needs no reference.  The grounded fraction is set to exactly 0 and
      ! exactly 1 away from the contour, so the test needs no tolerance.
      agr_eta = agree_det
      if (hybrid) then
        ! Try the tilt Laplacian first and keep the cell unless it FAILS.  Failure is
        ! either that it finds no stencil at all - it reads three cells, but the reference
        ! it is held against reads five - or that the reach exemption refuses to trust the
        ! reference it did build, which is what happens where the window straddles the
        ! flotation contour.  Both of those leave the cell undamped.  Stencil agreement
        ! needs no reference and no run of five, so it takes exactly those cells, and the
        ! stronger detector keeps every other one.
        call dg_tilt_detector_1d(tv, bv, hv, okv, gv, CS%dg_damp_kink_ref, &
                                 CS%dg_damp_kink_fit_tol, &
                                 A_h, A_bed, A_ref, href_d, det_ok, dmode, kink_c)
        agr_eta = .true.
        if (det_ok) then
          gwt = kink_c + (1.0 - kink_c) * &
                dg_gl_reach_wt(gv, max(dmode,1), CS%dg_damp_gl_reach)
          agr_eta = (gwt <= 0.0)
        endif
      endif
      if (agr_eta) then
        call dg_agree_1d(tv, bv, gv, hv, ilv, dy_cell, okv, CS%dg_damp_single_rule, &
                         CS%dg_damp_edge_rule, CS%dg_damp_reduce, CS%dg_damp_fields, &
                         A_dmp, A_spread, A_cen, n_cen, det_ok)
        if (CS%dg_damp_centred_amp .and. (A_dmp /= 0.0) .and. (n_cen > 0)) A_dmp = A_cen
        excess = abs(A_dmp) ; gwt = 1.0 ; href_d = href
      else
        ! Under the hybrid the detector ran above and gwt is already in hand; its answer
        ! is what chose this branch.
        if (.not.hybrid) then
          call dg_tilt_detector_1d(tv, bv, hv, okv, gv, CS%dg_damp_kink_ref, &
                                   CS%dg_damp_kink_fit_tol, &
                                   A_h, A_bed, A_ref, href_d, det_ok, dmode, kink_c)
          gwt = kink_c + (1.0 - kink_c) * &
                dg_gl_reach_wt(gv, max(dmode,1), CS%dg_damp_gl_reach)
        endif
        A_bed = f_gnd(i,j)*A_bed
        if (CS%dg_damp_excess_only) then
          A_dmp = dg_unexplained(A_h, A_ref, A_bed)
          excess = abs(A_dmp)
        else
          floor_A = max(abs(A_bed), abs(A_ref))
          excess = abs(A_h) - floor_A
          A_dmp = A_h
        endif
      endif
      if (det_ok .and. (excess > 0.0)) then
        dsurf = 0.0
        if (slope_gate) dsurf = 0.5*(sbar_c(i,j+1) - sbar_c(i,j-1))
        gate_now = CS%dg_damp_gate_form
        if (hybrid .and. agr_eta) gate_now = DAMP_GATE_AGREEMENT
        gam = gwt * dg_damp_gate_frac(gate_now, dsdh, excess, A_spread, href_d, &
                                      dsurf, CS%dg_tilt_damp_r_hi, CS%dg_damp_rho_s, &
                                      CS%dg_damp_slope_floor*dy_cell, CS%dg_damp_rho_g)
        kwant = gam * itau_eta
        kap = min(max(kwant - knat_eta, 0.0), rate_cap)
        r_eta(i,j) = kap*A_dmp
        if (diag_on) then
          d_ahat = d_ahat + (A_dmp**2) ; d_spread = max(d_spread, A_spread)
          d_gate = max(d_gate, gam) ; d_want = d_want + kwant
          d_got = d_got + kap
        endif
      endif
    endif

    ! --- xy-twist ---
    if (CS%dg_twist_damp .and. (.not.skip_w)) then
      ! Agreement reads the cross only, five cells along each axis, so the other 72 of the 81
      ! need not be gathered.  Today's twist reference uses the full square.
      ! Neither branch needs the stencils cleared first: today's writes all 81 entries, and
      ! agreement reads the cross alone, which the sweep below writes.  The two halves of that
      ! sweep are the image of each other under a quarter turn.
      if (agree_det) then
        do k = kglo, kghi
          kk = min(max(i+k, isd), ied) ; jj = min(max(j+k, jsd), jed)
          okw(k,0) = ice_ok(kk,j) .and. ((i+k >= isd) .and. (i+k <= ied))
          ww(k,0) = w_c(kk,j) ; bw(k,0) = b_w(kk,j) ; fw(k,0) = f_gnd(kk,j)
          gwx(k) = f_gnd(kk,j) ; iwx(k) = G%IdxT(kk,j)
          okw(0,k) = ice_ok(i,jj) .and. ((j+k >= jsd) .and. (j+k <= jed))
          ww(0,k) = w_c(i,jj) ; bw(0,k) = b_w(i,jj) ; fw(0,k) = f_gnd(i,jj)
          gwy(k) = f_gnd(i,jj) ; iwy(k) = G%IdyT(i,jj)
        enddo
      else
        do kj = kglo, kghi ; do k = kglo, kghi
          kk = min(max(i+k, isd), ied) ; jj = min(max(j+kj, jsd), jed)
          okw(k,kj) = ice_ok(kk,jj) .and. ((i+k >= isd) .and. (i+k <= ied)) &
                                    .and. ((j+kj >= jsd) .and. (j+kj <= jed))
          ww(k,kj) = w_c(kk,jj) ; bw(k,kj) = b_w(kk,jj)
          fw(k,kj) = f_gnd(kk,jj) ; hw(k,kj) = hbar_c(kk,jj)
          if (kj == 0) then
            gwx(k) = f_gnd(kk,jj) ; iwx(k) = G%IdxT(kk,j)
          endif
          if (k == 0) then
            gwy(kj) = f_gnd(kk,jj) ; iwy(kj) = G%IdyT(i,jj)
          endif
        enddo ; enddo
      endif

      ! The twist's reference reads a two-dimensional window, so the whole window has to
      ! stay on one side of the flotation contour before the tilt Laplacian may have it.
      agr_w = agree_det
      if (hybrid) then
        ngr = 0 ; nfl = 0
        do kj_g = -2, 2 ; do ki = -2, 2
          if (okw(ki,kj_g)) then
            if (fw(ki,kj_g) > 0.0) ngr = ngr + 1
            if (fw(ki,kj_g) < 1.0) nfl = nfl + 1
          endif
        enddo ; enddo
        agr_w = .not.((ngr == 0) .or. (nfl == 0))
      endif
      if (.not.agr_w) then

      ! A_w = (A_xi(w) + A_eta(w))/2, each axis with its own bias. All admissible bias pairs in
      ! the least-biased tier are averaged, so the choice does not depend on trial order.
      ntw = 0
      do tier = 1, 3
        do m = tw_lo(tier), tw_hi(tier)
          bx = bx_try(m) ; by = by_try(m)
          if (all(okw(Slo(bx):Shi(bx), Flo(by):Fhi(by))) .and. &
              all(okw(Flo(bx):Fhi(bx), Slo(by):Shi(by)))) then
            ntw = ntw + 1 ; bx_use(ntw) = bx ; by_use(ntw) = by
          endif
        enddo
        if (ntw > 0) exit
      enddo
      tw_ok = (ntw > 0)

      excess = -1.0
      if (tw_ok) then
        ! Blend the reference toward the kink's twist by the grounded-fraction span.
        kw_ok = .false. ; ksig = 0.0 ; kconf = 0.0
        if (CS%dg_damp_kink_ref) &
          call dg_kink_plane_2d(hw, fw, okw, -2, 2, CS%dg_damp_kink_tol, &
                                kqx, kqy, kq0, ksig, kconf, kw_ok)
        ksig = ksig * kconf
        A_w = 0.0 ; A_w_bed = 0.0 ; A_w_ref = 0.0 ; href_w = 0.0 ; reach_w = 1.0
        do itw = 1, ntw
          bx = bx_use(itw) ; by = by_use(itw)
          wrx(:) = 0.0 ; wry(:) = 0.0
          do k = -2, 2
            wrx(k) = dg_wref_at(hw, k, 0, bx, by)
            wry(k) = dg_wref_at(hw, 0, k, bx, by)
            if (kw_ok) then
              wrx(k) = ((1.0-ksig)*wrx(k)) + (ksig*dg_kink_twist_at(kqx, kqy, kq0, k, 0))
              wry(k) = ((1.0-ksig)*wry(k)) + (ksig*dg_kink_twist_at(kqx, kqy, kq0, 0, k))
            endif
          enddo
          ! Response 1 - (cos(theta_x) + cos(theta_y))/2: zero on a uniform twist,
          ! 2 on the checkerboard that no face jump can see.
          A_w     = A_w     + 0.5*(dg_d2_biased(ww(:,0), bx) + dg_d2_biased(ww(0,:), by))
          A_w_bed = A_w_bed + 0.5*(dg_d2_biased(bw(:,0), bx) + dg_d2_biased(bw(0,:), by))
          A_w_ref = A_w_ref + 0.5*(dg_d2_biased(wrx, bx) + dg_d2_biased(wry, by))
          href_w  = href_w  + (hw(0,0) + ((hw(Elo(bx),0) + hw(Ehi(bx),0)) + &
                                          (hw(0,Elo(by)) + hw(0,Ehi(by))))) / 5.0
          ! The reach exemption is a protection, so take the worst over the pairs
          ! rather than their mean.
          reach_w = min(reach_w, min(dg_gl_reach_wt(gwx, bx, CS%dg_damp_gl_reach), &
                                     dg_gl_reach_wt(gwy, by, CS%dg_damp_gl_reach)))
        enddo
        ! Averaged over the pairs TOGETHER, so that a smooth solution, on which
        ! every pair gives A_w == A_w_ref, still cancels term by term.
        inv_ntw = 1.0 / real(ntw)
        A_w = A_w * inv_ntw ; A_w_ref = A_w_ref * inv_ntw
        A_w_bed = f_gnd(i,j) * (A_w_bed * inv_ntw)
        href_w = href_w * inv_ntw
        if (CS%dg_damp_excess_only) then
          A_dmp = dg_unexplained(A_w, A_w_ref, A_w_bed)
          excess = abs(A_dmp)
        else
          excess = abs(A_w) - max(abs(A_w_bed), abs(A_w_ref))
          A_dmp = A_w
        endif
      endif
        ! The same failure test as the slopes: no admissible pair of stencils on both axes
        ! at once, or an exemption that refuses to trust the reference they built.
        if (hybrid .and. ((.not.tw_ok) .or. (reach_w <= 0.0))) agr_w = .true.
      endif
      if (agr_w) then
        ! Both axes, never a diagonal: the checkerboard alternates along x and along y, but is
        ! constant along either diagonal, so a diagonal stencil would read zero.  The second
        ! axis is also why the twist needs no single-stencil rule.
        call dg_agree_2d(ww, bw, fw, iwx, iwy, dx_cell, dy_cell, okw, CS%dg_damp_edge_rule, &
                         CS%dg_damp_reduce, CS%dg_damp_fields, &
                         A_dmp, A_spread, A_cen, n_cen, tw_ok)
        if (CS%dg_damp_centred_amp .and. (A_dmp /= 0.0) .and. (n_cen > 0)) A_dmp = A_cen
        excess = abs(A_dmp) ; gwt = 1.0 ; href_w = href
      endif

      if (excess > 0.0) then
        if (.not.agr_w) then
          ! Protection interpolated on the kink-fit confidence.
          ! The twist's second difference separates by axis, so the protection
          ! does too: each axis contributes the minimum over its own bias's reach.
          ! Where the reference carries the break there is nothing for the
          ! exemption to work around and the cell damps like any other; where it
          ! does not, the exemption is all there is.  Interpolated by confidence,
          ! not switched on it: a hard test at confidence zero would hand a cell
          ! FULL damping the moment the fit became marginally credible, while the
          ! blend had already reverted the reference to the uncorrected one --
          ! full strength against the worst floor, which is the combination this
          ! whole exercise exists to avoid.
          gwt = kconf + ((1.0 - kconf) * reach_w)
        endif
        dsurf = 0.0
        if (slope_gate) &
          dsurf = 0.25*((sbar_c(i+1,j+1) + sbar_c(i-1,j-1)) - &
                        (sbar_c(i-1,j+1) + sbar_c(i+1,j-1)))
        gate_now = CS%dg_damp_gate_form
        if (hybrid .and. agr_w) gate_now = DAMP_GATE_AGREEMENT
        gam = gwt * dg_damp_gate_frac(gate_now, dsdh, excess, A_spread, href_w, &
                                      dsurf, CS%dg_tilt_damp_r_hi, CS%dg_damp_rho_s, &
                                      CS%dg_damp_slope_floor*sqrt(dx_cell*dy_cell), &
                                      CS%dg_damp_rho_g)
        ! Credit the slower of the two sweeps.
        ! The checkerboard twist alternates along BOTH axes, so it survives as
        ! long as either sweep is slow: credit the transport with the smaller of
        ! the two rates, not their sum.
        kwant = gam * itau_w
        kap = min(max(kwant - min(knat_xi, knat_eta), 0.0), rate_cap)
        ! The corner pattern that carries this away is +,-,-,+ times a quarter of it, which
        ! changes w by -kap*A_dmp and leaves the tilts and the unweighted corner mean alone.
        r_w(i,j) = kap*A_dmp
        if (diag_on) then
          d_ahat = d_ahat + (A_dmp**2) ; d_spread = max(d_spread, A_spread)
          d_gate = max(d_gate, gam) ; d_want = d_want + kwant
          d_got = d_got + kap
        endif
      endif
    endif

    if (diag_on) then
      CS%dg_damp_ahat(i,j) = sqrt(d_ahat)
      CS%dg_damp_spread(i,j) = d_spread
      CS%dg_damp_gate(i,j) = d_gate
      CS%dg_damp_want(i,j) = d_want ; CS%dg_damp_got(i,j) = d_got
    endif
  enddo ; enddo

  ! Turn the removals into nodal rates.  The zigzag changes sign from each cell to the next, so
  ! the field taken away must do the same.  Its sign does, but its size does not: the detector
  ! and the rate both change from cell to cell, so the removal is an alternating pattern times a
  ! slowly changing envelope, and such a product also moves the smooth field.  The optional
  ! 1-2-1 high pass leaves a pure zigzag exactly as it is and takes out everything that changes
  ! linearly from cell to cell.  It is the same second difference the centred stencil uses,
  ! applied to the removal in place of the slopes.
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    rf_xi = r_xi(i,j) ; rf_eta = r_eta(i,j) ; rf_w = r_w(i,j)
    if (CS%dg_damp_filter) then
      ! A slope is filtered along its own direction, the twist along both, which is the same
      ! rule the stencils follow.  Where a neighbour is unusable the raw removal stands.
      if (ice_ok(i-1,j) .and. ice_ok(i+1,j)) &
        rf_xi = 0.25*((2.0*r_xi(i,j)) - (r_xi(i-1,j) + r_xi(i+1,j)))
      if (ice_ok(i,j-1) .and. ice_ok(i,j+1)) &
        rf_eta = 0.25*((2.0*r_eta(i,j)) - (r_eta(i,j-1) + r_eta(i,j+1)))
      if (all(ice_ok(i-1:i+1,j-1:j+1))) &
        rf_w = 0.0625*(((4.0*r_w(i,j)) - &
                        (2.0*((r_w(i-1,j) + r_w(i+1,j)) + (r_w(i,j-1) + r_w(i,j+1))))) + &
                       ((r_w(i-1,j-1) + r_w(i+1,j+1)) + (r_w(i-1,j+1) + r_w(i+1,j-1))))
    endif

    T_node(i,j,1,1) = ((0.5*rf_xi) + (0.5*rf_eta)) - (0.25*rf_w)
    T_node(i,j,1,2) = ((0.5*rf_xi) - (0.5*rf_eta)) + (0.25*rf_w)
    T_node(i,j,2,1) = ((0.5*rf_eta) - (0.5*rf_xi)) + (0.25*rf_w)
    T_node(i,j,2,2) = (-(0.5*rf_xi) - (0.5*rf_eta)) - (0.25*rf_w)

    ! The damper changes only the tilts and the twist, so it must move no mass.  Each
    ! increment above adds one sign to two corners and the other sign to the other two,
    ! which has no cell mean only where the two face lengths across the cell agree.  A
    ! latitude-longitude or tripolar grid does not give that, because corner_cell_weights
    ! builds the weights from the four face lengths, so remove whatever mean survives.
    ! The four corners are written out rather than swept by a 2x2 loop, which the
    ! compiler vectorizes into fused products that a grid rotation would not reproduce.
    Tbar = nodal_cell_mean(T_node(i,j,:,:), CS%cell_mean_w(i,j,:,:))
    T_node(i,j,1,1) = T_node(i,j,1,1) - Tbar
    T_node(i,j,2,2) = T_node(i,j,2,2) - Tbar
    T_node(i,j,2,1) = T_node(i,j,2,1) - Tbar
    T_node(i,j,1,2) = T_node(i,j,1,2) - Tbar

    if (diag_on) &
      CS%dg_damp_tend(i,j) = sqrt((((rf_xi**2) + (rf_eta**2)) / 12.0) + ((rf_w**2) / 144.0))
  enddo ; enddo

end subroutine dg_nodal_mode_damp_rate

!> Nodal rate that removes the alternating part of the DG(1) tilts and twist, measured against the
!! tilts the neighbouring cell means imply, where the artificial viscosity acts.
!!
!! The viscosity's jump penalty drives the thickness toward a continuous field with the cell means
!! it has.  In one dimension the corner values of such a field obey n(k+1) = 2*hbar(k) - n(k), so any
!! misfit at one corner is carried from cell to cell with its sign flipped and no loss, and it shows
!! as a tilt that alternates between cells without a face jump.  Each tilt is compared with the
!! mean-supported one, e = t - t_ref, and only the part of e that alternates, e - (e_- + e_+)/2, is
!! removed.  A smooth difference between the two is left alone: along a flow line the upwind DG tilt
!! legitimately differs from a centred difference of the means, and relaxing that difference would
!! change the fluxes.  The reference is the slope of a line through the cell's own mean fitted by
!! least squares to its neighbours' means along the axis, which is the centred difference in the
!! interior and one-sided beside an ice edge.  It is read in surface form, so grounded ice that
!! follows its bed is left alone, and only from neighbours of the same flotation state, so it never
!! reaches across a grounding line.  A cell the grounding line crosses is neither relaxed nor used as
!! a neighbour: the real change of slope is there.  The exception is DG1_TILT_RELAX_FREE_EDGE, for
!! an axis with no artificial-viscosity face on one side: that node has no pin at all, so a bounded
!! contamination from reading across the grounding line is better than leaving it to drift.
!! The rate along each axis is a fraction of the
!! viscosity's own jump decay rate on the cell's two faces on that axis, so the removal acts where
!! and as fast as the penalty that feeds the alternating field.
subroutine dg1_tilt_relax_rate(CS, G, hmask, h_nodal_in, dt, R_node)
  type(ice_shelf_dyn_CS), intent(in)  :: CS   !< Ice shelf dynamics control structure
  type(ocean_grid_type),  intent(inout) :: G  !< Ocean grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), intent(in) :: hmask !< Cell mask
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(in)  :: h_nodal_in !< Corner thicknesses [Z ~> m]
  real,                   intent(in)  :: dt   !< Time step [T ~> s]
  real, dimension(SZDI_(G),SZDJ_(G),2,2), intent(out) :: R_node !< Nodal relaxation rate [Z T-1 ~> m s-1]

  integer, dimension(SZDI_(G),SZDJ_(G)) :: state ! 1 grounded, -1 floating, 0 crossed or not ice
  real, dimension(SZDI_(G),SZDJ_(G)) :: gx, gy ! Viscosity jump decay rate along each axis, the
                                               ! larger of the cell's two faces [T-1 ~> s-1]
  real, dimension(SZDI_(G),SZDJ_(G)) :: e_xi, e_eta, e_w ! Tilt and twist less their
                                               ! mean-supported values [Z ~> m]
  logical, dimension(SZDI_(G),SZDJ_(G)) :: ok_xi, ok_eta, ok_w ! Those residuals exist
  real :: gxs, gys   ! The rate on each axis for this cell [T-1 ~> s-1]
  real :: frac_use   ! Multiplier on the axis rate: 1 under the speed law [nondim]
  real :: unat_x, unat_y ! Credited part of the cell-centred speed on each axis [L T-1 ~> m s-1]
  real :: a_xi, a_eta, a_w ! Alternating part of each residual [Z ~> m]
  real :: r_xi, r_eta, r_w ! Removal rates per mode [Z T-1 ~> m s-1]
  real :: k_xi, k_eta, k_w ! Delivered rates per mode [T-1 ~> s-1]
  real :: num, den  ! Numerator [Z L ~> m2] and denominator [L2 ~> m2] of the least-squares slope
  real :: dm, dp    ! Distances to the neighbours behind and ahead along the axis [L ~> m]
  real :: esum      ! Sum of the usable neighbour residuals [Z ~> m]
  real :: esum_x, esum_y ! The same along each axis, for the twist [Z ~> m]
  integer :: n_x, n_y, nnb ! Numbers of usable neighbours
  real :: rate_cap  ! Largest rate the explicit step may carry [T-1 ~> s-1]
  real :: f         ! Grounded fraction clipped to [0,1] [nondim]
  real :: Tbar      ! Cell mean of the nodal increment, removed so no mass moves [Z T-1 ~> m s-1]
  logical :: use_m, use_p ! The neighbour behind or ahead is usable
  logical :: diag_on
  ! For the unpinned cells the grounding line crosses, whose tilts are read in surface elevation.
  real, dimension(SZDI_(G),SZDJ_(G)) :: sbar   ! Mean of the corner surface elevations [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)) :: es_xi, es_eta, es_w ! Surface tilts and twist less their
                                               ! mean-supported values [Z ~> m]
  logical, dimension(SZDI_(G),SZDJ_(G)) :: ice ! The cell holds ice and is inside the outer ring
  logical, dimension(SZDI_(G),SZDJ_(G)) :: oks_xi, oks_eta, oks_w ! Those surface residuals exist
  logical, dimension(SZDI_(G),SZDJ_(G)) :: fe_x, fe_y ! The cell has no artificial-viscosity face on
                                               ! one side of this axis, so that node has no pin
  logical :: wx, wy ! Any ice neighbour is usable on this axis, for this cell
  logical :: wide   ! A neighbour is usable whenever it holds ice, whatever its flotation state
  real :: ts_w      ! The cell's surface twist [Z ~> m]
  real :: m_cell    ! Thickness change per unit surface change for the cell [nondim]
  real :: s_c(2,2)  ! Surface elevation at the cell's corners [Z ~> m]
  real :: fg_cell   ! Grounded area fraction of a crossed cell, clipped to [0,1] [nondim]
  real :: dhds_cell ! Thickness change per unit surface change for the whole cell, the harmonic
                    ! mean of the corner values, which keeps the removal a pure tilt [nondim]
  real :: ts_xi, ts_eta ! The cell's surface tilts [Z ~> m]
  real :: rr        ! Ratio of ice to ocean density [nondim]
  logical :: pin_x, pin_y ! The cell has a relaxed neighbour on both sides along that axis
  integer :: m, n
  real, parameter :: f_eps = 1.0e-6 ! Grounded-fraction tolerance for a pure flotation state [nondim]
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  R_node(isc:iec,jsc:jec,:,:) = 0.0
  diag_on = (CS%id_dg_tilt_relax_rate > 0) .or. (CS%id_dg_tilt_relax_tend > 0) .or. &
            (CS%id_dg_tilt_relax_alt_x > 0) .or. (CS%id_dg_tilt_relax_alt_y > 0) .or. &
            (CS%id_dg_tilt_relax_alt_w > 0) .or. (CS%id_dg_tilt_relax_state > 0)
  if (diag_on) then
    CS%dg_tilt_relax_rate(:,:) = 0.0 ; CS%dg_tilt_relax_tend(:,:) = 0.0
    CS%dg_tilt_relax_alt_x(:,:) = 0.0 ; CS%dg_tilt_relax_alt_y(:,:) = 0.0
    CS%dg_tilt_relax_alt_w(:,:) = 0.0 ; CS%dg_tilt_relax_state(:,:) = 0.0
  endif
  if (.not.CS%dg_tilt_relax) return
  if (.not.associated(CS%bed_node)) call MOM_error(FATAL, &
    "dg1_tilt_relax_rate: DG1_TILT_RELAX requires a nodal bed (CS%bed_node).")
  ! The residuals are needed one cell beyond the compute domain and read means one further, and the
  ! outermost ring of the data domain is skipped.
  if ((G%isc - G%isd < 3) .or. (G%jsc - G%jsd < 3)) call MOM_error(FATAL, &
    "dg1_tilt_relax_rate: DG1_TILT_RELAX needs a halo of at least 3 cells; increase NIHALO/NJHALO.")

  rate_cap = 0.5 / max(dt, tiny(dt))

  ! Flotation state on the data domain.  The comparison variable itself is the corner surface,
  ! built where the residuals are; only the state is needed here, to decide which neighbours a
  ! cell may use.
  do j = jsd, jed ; do i = isd, ied
    state(i,j) = 0
    if (hmask(i,j) /= 1.0) cycle
    if ((i == isd) .or. (i == ied) .or. (j == jsd) .or. (j == jed)) cycle
    f = min(max(CS%ground_frac(i,j), 0.0), 1.0)
    if (f >= 1.0 - f_eps) then
      state(i,j) = 1
    elseif (f <= f_eps) then
      state(i,j) = -1
    endif
  enddo ; enddo

  ! An axis has a free edge where one side carries no artificial-viscosity face: dg1_art_visc_face
  ! is only built where both cells are hmask 1 or 3, so beyond a domain edge or a calving front
  ! there is none and the outer node has no pin at all.  A thickness boundary does NOT count: the
  ! face is active there and the node is held.
  do j = jsd+1, jed-1 ; do i = isd+1, ied-1
    fe_x(i,j) = .not.(((hmask(i-1,j) == 1.0) .or. (hmask(i-1,j) == 3.0)) .and. &
                      ((hmask(i+1,j) == 1.0) .or. (hmask(i+1,j) == 3.0)))
    fe_y(i,j) = .not.(((hmask(i,j-1) == 1.0) .or. (hmask(i,j-1) == 3.0)) .and. &
                      ((hmask(i,j+1) == 1.0) .or. (hmask(i,j+1) == 3.0)))
  enddo ; enddo
  fe_x(isd,:) = .false. ; fe_x(ied,:) = .false. ; fe_x(:,jsd) = .false. ; fe_x(:,jed) = .false.
  fe_y(isd,:) = .false. ; fe_y(ied,:) = .false. ; fe_y(:,jsd) = .false. ; fe_y(:,jed) = .false.

  ! A neighbour is usable whatever its flotation state under the unpinned rule, which reaches a
  ! crossed cell on purpose, and under the wide rule.  Otherwise it must share this cell's state,
  ! which is what makes the surface form reproduce the thickness form it replaces.
  wide = CS%dg_tilt_relax_wide_nb
  rr = CS%density_ice / CS%density_ocean_avg
  sbar(:,:) = 0.0 ; ice(:,:) = .false.
  do j = jsd+1, jed-1 ; do i = isd+1, ied-1
    if (hmask(i,j) /= 1.0) cycle
    ice(i,j) = .true.
    do n = 1, 2 ; do m = 1, 2
      s_c(m,n) = max(h_nodal_in(i,j,m,n) - CS%bed_node(i-2+m,j-2+n), (1.0-rr)*h_nodal_in(i,j,m,n))
    enddo ; enddo
    sbar(i,j) = 0.25*((s_c(1,1) + s_c(2,2)) + (s_c(2,1) + s_c(1,2)))
  enddo ; enddo
  es_xi(:,:) = 0.0 ; es_eta(:,:) = 0.0 ; es_w(:,:) = 0.0
  oks_xi(:,:) = .false. ; oks_eta(:,:) = .false. ; oks_w(:,:) = .false.
  do j = jsc-1, jec+1 ; do i = isc-1, iec+1
    if (.not.ice(i,j)) cycle
    ! An axis is widened to any ice neighbour under the wide rule, or where it has a free edge
    ! and the free-edge rule is on: there the only alternative is no reference at all.
    wx = wide .or. (CS%dg_tilt_relax_free_edge .and. fe_x(i,j))
    wy = wide .or. (CS%dg_tilt_relax_free_edge .and. fe_y(i,j))
    ! Without either, a crossed cell is left out entirely, exactly as the thickness form leaves
    ! it out, so that it is neither relaxed nor offered to a neighbour as a reference.
    if ((state(i,j) == 0) .and. (.not.(wx .or. wy))) cycle
    do n = 1, 2 ; do m = 1, 2
      s_c(m,n) = max(h_nodal_in(i,j,m,n) - CS%bed_node(i-2+m,j-2+n), (1.0-rr)*h_nodal_in(i,j,m,n))
    enddo ; enddo
    ts_xi  = 0.5*((s_c(2,1) - s_c(1,1)) + (s_c(2,2) - s_c(1,2)))
    ts_eta = 0.5*((s_c(1,2) - s_c(1,1)) + (s_c(2,2) - s_c(2,1)))
    use_m = ice(i-1,j) ; use_p = ice(i+1,j)
    if (.not.wx) then
      use_m = use_m .and. (state(i-1,j) == state(i,j))
      use_p = use_p .and. (state(i+1,j) == state(i,j))
    endif
    if (use_m .or. use_p) then
      dm = G%dxCu(I-1,j) ; dp = G%dxCu(I,j)
      if (use_m .and. use_p) then
        num = (dp*(sbar(i+1,j) - sbar(i,j))) + (dm*(sbar(i,j) - sbar(i-1,j)))
        den = (dp*dp) + (dm*dm)
      elseif (use_p) then
        num = dp*(sbar(i+1,j) - sbar(i,j)) ; den = dp*dp
      else
        num = dm*(sbar(i,j) - sbar(i-1,j)) ; den = dm*dm
      endif
      es_xi(i,j) = ts_xi - ((num / den) * G%dxT(i,j)) ; oks_xi(i,j) = .true.
    endif
    use_m = ice(i,j-1) ; use_p = ice(i,j+1)
    if (.not.wy) then
      use_m = use_m .and. (state(i,j-1) == state(i,j))
      use_p = use_p .and. (state(i,j+1) == state(i,j))
    endif
    if (use_m .or. use_p) then
      dm = G%dyCv(i,J-1) ; dp = G%dyCv(i,J)
      if (use_m .and. use_p) then
        num = (dp*(sbar(i,j+1) - sbar(i,j))) + (dm*(sbar(i,j) - sbar(i,j-1)))
        den = (dp*dp) + (dm*dm)
      elseif (use_p) then
        num = dp*(sbar(i,j+1) - sbar(i,j)) ; den = dp*dp
      else
        num = dm*(sbar(i,j) - sbar(i,j-1)) ; den = dm*dm
      endif
      es_eta(i,j) = ts_eta - ((num / den) * G%dyT(i,j)) ; oks_eta(i,j) = .true.
    endif

    ! The twist, from the cross difference of the four diagonal neighbours' mean surfaces.  It
    ! reads means rather than neighbour tilts, so a tilt mode cannot leak into its reference.
    if (CS%dg_tilt_relax_twist) then
      use_m = (ice(i-1,j-1) .and. ice(i+1,j+1)) .and. (ice(i-1,j+1) .and. ice(i+1,j-1))
      use_m = use_m .and. (state(i,j) /= 0)
      if (.not.wide) use_m = use_m .and. &
        (((state(i-1,j-1) == state(i,j)) .and. (state(i+1,j+1) == state(i,j))) .and. &
         ((state(i-1,j+1) == state(i,j)) .and. (state(i+1,j-1) == state(i,j))))
      if (use_m) then
        ts_w = (s_c(2,2) - s_c(1,2)) - (s_c(2,1) - s_c(1,1))
        es_w(i,j) = ts_w - (0.25*((sbar(i+1,j+1) + sbar(i-1,j-1)) - &
                                  (sbar(i-1,j+1) + sbar(i+1,j-1))))
        oks_w(i,j) = .true.
      endif
    endif
  enddo ; enddo

  ! The surface residuals are the residuals.  The bed tilts are already inside them, because
  ! s_c subtracts the bed at each grounded corner, so no separate bed stencil is needed.
  e_xi(:,:) = es_xi(:,:) ; ok_xi(:,:) = oks_xi(:,:)
  e_eta(:,:) = es_eta(:,:) ; ok_eta(:,:) = oks_eta(:,:)
  e_w(:,:) = es_w(:,:) ; ok_w(:,:) = oks_w(:,:)

  ! Rate along each axis.  With DG1_TILT_RELAX_U_CUT set this is the speed law
  ! max(u_cut - u_credit*|u_n|, 0)/dx, which is grid-invariant and carries no dependence on the
  ! artificial viscosity: that gate reads the face jump, which this mode does not produce, so it is
  ! smallest exactly where the error is purest.  Otherwise it is the legacy fraction of the larger
  ! viscosity jump decay rate of the cell's two faces on that axis, 4*amp*nu/dx_perp^2, from the
  ! face viscosities the spatial operator wrote on this same stage.  The cell-centred speed sums its
  ! four corners in diagonal pairs, so a quarter turn of the grid permutes the same additions.
  gx(:,:) = 0.0 ; gy(:,:) = 0.0
  frac_use = CS%dg_tilt_relax_frac
  if (CS%dg_tilt_relax_u_cut > 0.0) frac_use = 1.0
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    if (CS%dg_tilt_relax_u_cut > 0.0) then
      unat_x = CS%dg_tilt_relax_u_credit * 0.25 * &
               ((abs(CS%u_shelf(I-1,J-1)) + abs(CS%u_shelf(I,J))) + &
                (abs(CS%u_shelf(I,J-1))   + abs(CS%u_shelf(I-1,J))))
      unat_y = CS%dg_tilt_relax_u_credit * 0.25 * &
               ((abs(CS%v_shelf(I-1,J-1)) + abs(CS%v_shelf(I,J))) + &
                (abs(CS%v_shelf(I,J-1))   + abs(CS%v_shelf(I-1,J))))
      gx(i,j) = max(CS%dg_tilt_relax_u_cut - unat_x, 0.0) * G%IdxT(i,j)
      gy(i,j) = max(CS%dg_tilt_relax_u_cut - unat_y, 0.0) * G%IdyT(i,j)
    else
      gx(i,j) = (4.0*DG1_WB_JUMP_RATE_AMP) * &
                max(CS%dg_art_visc_nu_u(I-1,j) / (G%dxCu(I-1,j)*G%dxCu(I-1,j)), &
                    CS%dg_art_visc_nu_u(I,j) / (G%dxCu(I,j)*G%dxCu(I,j)))
      gy(i,j) = (4.0*DG1_WB_JUMP_RATE_AMP) * &
                max(CS%dg_art_visc_nu_v(i,J-1) / (G%dyCv(i,J-1)*G%dyCv(i,J-1)), &
                    CS%dg_art_visc_nu_v(i,J) / (G%dyCv(i,J)*G%dyCv(i,J)))
    endif
  enddo ; enddo

  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    if (diag_on) CS%dg_tilt_relax_state(i,j) = real(state(i,j))
    gxs = gx(i,j) ; gys = gy(i,j)
    if ((gxs <= 0.0) .and. (gys <= 0.0)) cycle

    ! Which axes may act.  A pure cell acts on both.  A cell the grounding line crosses acts only
    ! where the wide rule is on, or where the axis has a free edge and nothing else holds its outer
    ! node.  Everything below is then common to every cell.
    wx = (state(i,j) /= 0) .or. wide .or. (CS%dg_tilt_relax_free_edge .and. fe_x(i,j))
    wy = (state(i,j) /= 0) .or. wide .or. (CS%dg_tilt_relax_free_edge .and. fe_y(i,j))
    if (.not.(wx .or. wy)) cycle

    r_xi = 0.0 ; r_eta = 0.0 ; r_w = 0.0 ; k_xi = 0.0 ; k_eta = 0.0 ; k_w = 0.0
    a_xi = 0.0 ; a_eta = 0.0 ; a_w = 0.0

    ! The alternating part of each residual: the residual less the mean of its usable neighbours
    ! along the same axis.  A cell with no usable neighbour cannot tell a zigzag from structure.
    if (((gxs > 0.0) .and. ok_xi(i,j)) .and. wx) then
      esum = 0.0 ; nnb = 0
      if (ok_xi(i-1,j) .and. ok_xi(i+1,j)) then
        esum = e_xi(i-1,j) + e_xi(i+1,j) ; nnb = 2
      elseif (ok_xi(i-1,j)) then
        esum = e_xi(i-1,j) ; nnb = 1
      elseif (ok_xi(i+1,j)) then
        esum = e_xi(i+1,j) ; nnb = 1
      endif
      if (nnb > 0) then
        a_xi = e_xi(i,j) - (esum / real(nnb))
        k_xi = min(frac_use*gxs, rate_cap)
        r_xi = k_xi * a_xi
      endif
    endif

    if (((gys > 0.0) .and. ok_eta(i,j)) .and. wy) then
      esum = 0.0 ; nnb = 0
      if (ok_eta(i,j-1) .and. ok_eta(i,j+1)) then
        esum = e_eta(i,j-1) + e_eta(i,j+1) ; nnb = 2
      elseif (ok_eta(i,j-1)) then
        esum = e_eta(i,j-1) ; nnb = 1
      elseif (ok_eta(i,j+1)) then
        esum = e_eta(i,j+1) ; nnb = 1
      endif
      if (nnb > 0) then
        a_eta = e_eta(i,j) - (esum / real(nnb))
        k_eta = min(frac_use*gys, rate_cap)
        r_eta = k_eta * a_eta
      endif
    endif

    ! The twist alternates along both axes, so its neighbours are the four sharing a face.  The two
    ! axis sums are formed first and then added, which a quarter turn only exchanges.
    if (ok_w(i,j)) then
      esum_x = 0.0 ; n_x = 0 ; esum_y = 0.0 ; n_y = 0
      if (ok_w(i-1,j) .and. ok_w(i+1,j)) then
        esum_x = e_w(i-1,j) + e_w(i+1,j) ; n_x = 2
      elseif (ok_w(i-1,j)) then
        esum_x = e_w(i-1,j) ; n_x = 1
      elseif (ok_w(i+1,j)) then
        esum_x = e_w(i+1,j) ; n_x = 1
      endif
      if (ok_w(i,j-1) .and. ok_w(i,j+1)) then
        esum_y = e_w(i,j-1) + e_w(i,j+1) ; n_y = 2
      elseif (ok_w(i,j-1)) then
        esum_y = e_w(i,j-1) ; n_y = 1
      elseif (ok_w(i,j+1)) then
        esum_y = e_w(i,j+1) ; n_y = 1
      endif
      if (n_x + n_y > 0) then
        a_w = e_w(i,j) - ((esum_x + esum_y) / real(n_x + n_y))
        k_w = min(frac_use*max(gxs, gys), rate_cap)
        r_w = k_w * a_w
      endif
    endif

    ! The residual was read in surface elevation, so convert the amount to thickness with ONE
    ! factor for the cell: the harmonic mean of dh/ds over it, which is 1 where the cell is wholly
    ! grounded and 1/(1-rr) where it is wholly afloat.  A factor that varied between the corners
    ! would turn the removal from a pure tilt into a tilt plus a twist, and the mode is a pure
    ! tilt of the THICKNESS nodes.  The same conversion dg1_wb_slope_mean makes across a face.
    ! The grounded AREA fraction keeps the factor continuous as the grounding line sweeps past a
    ! corner; "SEP2" gives a true area fraction, "SEP3" a counted one that moves in steps.
    fg_cell = min(max(CS%ground_frac(i,j), 0.0), 1.0)
    m_cell = 1.0 / (fg_cell + ((1.0 - fg_cell) * (1.0 - rr)))
    r_xi = r_xi * m_cell ; r_eta = r_eta * m_cell ; r_w = r_w * m_cell
    a_xi = a_xi * m_cell ; a_eta = a_eta * m_cell ; a_w = a_w * m_cell

    ! The corner pattern of the removal, as in dg_nodal_mode_damp_rate, with the cell mean taken out
    ! so that no mass moves on a grid whose opposite faces differ in length.
    R_node(i,j,1,1) = ((0.5*r_xi) + (0.5*r_eta)) - (0.25*r_w)
    R_node(i,j,1,2) = ((0.5*r_xi) - (0.5*r_eta)) + (0.25*r_w)
    R_node(i,j,2,1) = ((0.5*r_eta) - (0.5*r_xi)) + (0.25*r_w)
    R_node(i,j,2,2) = (-(0.5*r_xi) - (0.5*r_eta)) - (0.25*r_w)
    Tbar = nodal_cell_mean(R_node(i,j,:,:), CS%cell_mean_w(i,j,:,:))
    R_node(i,j,1,1) = R_node(i,j,1,1) - Tbar
    R_node(i,j,2,2) = R_node(i,j,2,2) - Tbar
    R_node(i,j,2,1) = R_node(i,j,2,1) - Tbar
    R_node(i,j,1,2) = R_node(i,j,1,2) - Tbar

    if (diag_on) then
      CS%dg_tilt_relax_rate(i,j) = max(max(k_xi, k_eta), k_w)
      CS%dg_tilt_relax_tend(i,j) = sqrt((((r_xi**2) + (r_eta**2)) / 12.0) + ((r_w**2) / 144.0))
      CS%dg_tilt_relax_alt_x(i,j) = a_xi ; CS%dg_tilt_relax_alt_y(i,j) = a_eta
      CS%dg_tilt_relax_alt_w(i,j) = a_w
    endif
  enddo ; enddo

end subroutine dg1_tilt_relax_rate


!> Advance CS%h_nodal one step with SSP-RK2, including the sources and mode damping, and return
!! the stage-averaged face fluxes.
subroutine ice_shelf_advect_DG1_nodal(CS, ISS, G, time_step, hmask, uh_ice, vh_ice)
  type(ice_shelf_dyn_CS), intent(inout) :: CS
  type(ice_shelf_state),  intent(in)    :: ISS
  type(ocean_grid_type),  intent(inout) :: G
  real,                   intent(in)    :: time_step
  real, dimension(SZDI_(G),SZDJ_(G)),       intent(inout) :: hmask
  real, dimension(SZDIB_(G),SZDJ_(G)),      intent(inout) :: uh_ice
  real, dimension(SZDI_(G),SZDJB_(G)),      intent(inout) :: vh_ice

  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: h0, h_curr, rhs
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: S_node ! Nodal source rate [Z T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: T_node ! Nodal mode-damping rate [Z T-1 ~> m s-1]
  real, dimension(SZDI_(G),SZDJ_(G),2,2) :: R_node ! Nodal tilt-relaxation rate [Z T-1 ~> m s-1]
  real, dimension(2,2) :: dh
  real :: tau_chk ! Shortest mode-damper relaxation time in the domain [T ~> s]
  real :: u_chk   ! Fastest removal speed either tilt term asks for in a cell [L T-1 ~> m s-1]
  character(len=24) :: n1, n2 ! Numbers for the warnings; the text is concatenated separately
  integer :: i, j, isc, iec, jsc, jec, isd, ied, jsd, jed, a, b

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  uh_ice(:,:) = 0.0 ; vh_ice(:,:) = 0.0

  ! Dirichlet thickness BC.  A boundary cell keeps the nodal thicknesses it was given, so that it
  ! can carry a profile along the inflow face; the stage update below already steps only hmask == 1
  ! cells, and recompute_h_shelf_from_nodal covers hmask == 3, so h_shelf follows those corners.
  ! Nothing is flattened here.  The former version snapped all four corners of every cell with a
  ! non-zero h_bdry_val to that one value, which both used a different test from the flux, the
  ! viscosity and the front term -- all of which key on hmask == 3 -- and left a boundary that was
  ! piecewise constant, so an interior neighbour could not hold a gradient against it.
  h0(:,:,:,:) = CS%h_nodal(:,:,:,:)
  call pass_corner_field(h0, G)

  ! A nodal source S enters as M*S, so after M^-1 it is added to dh/dt directly.
  call project_h_source_rate_to_nodes(CS, ISS, G, S_node)

  if (CS%debug) then
    call check_nodal_source_conservation(CS, ISS, G, CS%h_source_rate, S_node, "basal+surface")
    if (CS%dg_basal_source_sem2) call check_xi_basal_consistency(CS, ISS, G)
  endif

  ! Stage 1: positivity limit, spatial operator, M^-1, Euler step.
  call nodal_positivity_limit(CS, G, ISS)
  call DG1_nodal_spatial_operator(CS, G, hmask, CS%h_nodal, rhs, uh_ice, vh_ice, time_step)

  ! One floor on the delivered relaxation time, warned once. It is not a stability
  ! bound: the 0.5/dt rate cap covers stability outright. Once tau falls below the
  ! advection step that cap binds, and the DELIVERED relaxation time is the advection
  ! step rather than the tau that was asked for. Lowering tau further then changes
  ! nothing, so a tau ladder goes flat for a reason that has nothing to do with the
  ! solution. Sharp, and with no free constant in it.
  !
  ! A second floor, against the lagged velocity feedback, was tried here and removed.
  ! A steady state can carry a small grounding-line oscillation, because the velocity
  ! is frozen for one ICE_VELOCITY_TIMESTEP while the thickness keeps advancing.
  ! Measured on MISMIP+ at 10 km it is about 0.007 cells rms, some 70 m, with a
  ! 40-50 yr period, and it falls to 0.0003 cells when the velocity step is quartered.
  ! The tilt terms raise the loop gain and so bring it into view sooner, but they do
  ! not cause it: the same period appears with every tilt term off, and it is
  ! unaffected by the number of nonlinear substeps per velocity solve. No honest test
  ! could be written for it at this point. The clock that governs it is the velocity
  ! step over the cell traversal time dx/u, which had to stay under about 0.1 in that
  ! configuration -- but that number is a loop gain, so it moves with geometry,
  ! friction and buttressing; and this code runs once, on the first step, before the
  ! velocity field that sets the clock has developed. Any fixed threshold tested here
  ! would be arbitrary AND evaluated at the wrong time. If a fluctuating steady state
  ! does appear, shorten ICE_VELOCITY_TIMESTEP; lengthening the relaxation does not
  ! help, and was measured to make it worse.

  ! Warn once if the shortest relaxation time is below the advection step, where the rate cap
  ! binds. The tilt relaxation's speed law shares the clock, and at rest delivers U_CUT/dx.
  if ((CS%dg_tilt_damp .or. CS%dg_twist_damp .or. &
       (CS%dg_tilt_relax .and. (CS%dg_tilt_relax_u_cut > 0.0))) .and. &
      .not.CS%dg_tilt_damp_dt_warned) then
    ! tau = dx/(2*c*U_CUT) over the shortest cell dimension.
    tau_chk = huge(1.0)
    do j = jsc, jec ; do i = isc, iec
      u_chk = 0.0
      if (CS%dg_tilt_damp .or. CS%dg_twist_damp) u_chk = CS%dg_damp_advective_c * CS%dg_damp_u_cut
      if (CS%dg_tilt_relax) u_chk = max(u_chk, CS%dg_tilt_relax_u_cut)
      tau_chk = min(tau_chk, 0.5 / max(u_chk * max(G%IdxT(i,j), G%IdyT(i,j)), tiny(1.0)))
    enddo ; enddo
    ! Every PE reaches this collective.
    call min_across_PEs(tau_chk)
    if (is_root_pe()) then
      if (tau_chk < time_step) then
        write(n1,'(ES10.3)') tau_chk ; write(n2,'(ES10.3)') time_step
        call MOM_error(WARNING, "DG(1) tilt terms: shortest relaxation time "//&
             trim(adjustl(n1))//" s is below the advection step "//trim(adjustl(n2))//&
             " s. The rate cap binds, so the delivered relaxation time IS the advection "//&
             "step and not the one requested; lowering it further will change nothing.")
      endif
    endif
    CS%dg_tilt_damp_dt_warned = .true.
  endif

  call dg_nodal_mode_damp_rate(CS, G, hmask, CS%h_nodal, time_step, T_node)
  if (CS%dg_tilt_relax) then
    call dg1_tilt_relax_rate(CS, G, hmask, CS%h_nodal, time_step, R_node)
    do j = jsc, jec ; do i = isc, iec ; do b = 1, 2 ; do a = 1, 2
      T_node(i,j,a,b) = T_node(i,j,a,b) + R_node(i,j,a,b)
    enddo ; enddo ; enddo ; enddo
  endif
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    call apply_nodal_DG_mass_inverse(CS%Minv_nodal(i,j,:,:,:,:), &
                                     rhs(i,j,:,:), dh)
    do b = 1, 2 ; do a = 1, 2
      CS%h_nodal(i,j,a,b) = h0(i,j,a,b) &
                           + time_step * (dh(a,b) + S_node(i,j,a,b) + T_node(i,j,a,b))
    enddo ; enddo
  enddo ; enddo
  call pass_corner_field(CS%h_nodal, G)

  ! Stage 2: positivity limit, spatial operator, M^-1, SSP-RK2 combination.
  call nodal_positivity_limit(CS, G, ISS)
  h_curr(:,:,:,:) = CS%h_nodal(:,:,:,:)
  call DG1_nodal_spatial_operator(CS, G, hmask, h_curr, rhs, uh_ice, vh_ice, time_step)
  call dg_nodal_mode_damp_rate(CS, G, hmask, h_curr, time_step, T_node)
  if (CS%dg_tilt_relax) then
    call dg1_tilt_relax_rate(CS, G, hmask, h_curr, time_step, R_node)
    do j = jsc, jec ; do i = isc, iec ; do b = 1, 2 ; do a = 1, 2
      T_node(i,j,a,b) = T_node(i,j,a,b) + R_node(i,j,a,b)
    enddo ; enddo ; enddo ; enddo
  endif
  do j = jsc, jec ; do i = isc, iec
    if (hmask(i,j) /= 1.0) cycle
    call apply_nodal_DG_mass_inverse(CS%Minv_nodal(i,j,:,:,:,:), &
                                     rhs(i,j,:,:), dh)
    do b = 1, 2 ; do a = 1, 2
      CS%h_nodal(i,j,a,b) = 0.5*h0(i,j,a,b) &
                           + 0.5*(h_curr(i,j,a,b) &
                                  + time_step*(dh(a,b) + S_node(i,j,a,b) + T_node(i,j,a,b)))
    enddo ; enddo
  enddo ; enddo
  call pass_corner_field(CS%h_nodal, G)

  ! Keep the consumed source for diagnostics and reset the buffers.
  CS%h_source_rate_last(:,:) = CS%h_source_rate(:,:)
  CS%h_source_rate(:,:) = 0.0
  CS%h_source_rate_bmb(:,:) = 0.0

  ! Final positivity limit.
  call nodal_positivity_limit(CS, G, ISS)

  ! Average uh_ice, vh_ice over the 2 stages (SSP-RK2 equal weight).
  uh_ice(:,:) = 0.5 * uh_ice(:,:)
  vh_ice(:,:) = 0.5 * vh_ice(:,:)
  call pass_vector(uh_ice, vh_ice, G%domain, TO_ALL, CGRID_NE)
end subroutine ice_shelf_advect_DG1_nodal

end module MOM_ice_shelf_dynamics
