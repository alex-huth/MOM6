! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Implements the simplified Material Point Method (sMPM) for ice shelf dynamics,
!! following Huth et al. (2021) and the GIMPM reference implementation.
!!
!! Particles carry ice thickness, deformation history, and material properties
!! Lagrangianly, replacing the Eulerian advection of ISS%h_shelf.  The existing
!! MOM6 SSA (CG) solver is reused; MPM contributes by:
!!   1. Maintaining ISS%h_shelf and ISS%hmask via P2G (particles-to-grid)
!!   2. Updating particle state (F, H, position) via G2P after each velocity solve
!!   3. Splitting and migrating particles as the shelf front advances
module MOM_ice_shelf_MPM

use MOM_coms,          only : PE_here, num_PEs
use MOM_domains,       only : pass_var, CORNER
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, is_root_pe
use MOM_file_parser,   only : get_param, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_ice_shelf_state, only : ice_shelf_state

implicit none ; private

public :: MPM_CS, MPM_init, MPM_end
public :: MPM_P2G, MPM_G2P, MPM_Lagrangian_update
public :: MPM_split, MPM_migrate, MPM_enforce_1d_test
public :: MPM_update_masks
public :: MPM_compute_visc, MPM_compute_is_front_cell
public :: MPM_store_pre_solve_vel
public :: MPM_P2G_velocity, MPM_compute_basal_trac
public :: MPM_write_vtu, MPM_save_restart, MPM_restore_restart
public :: smpm_shape, smpm_grad, gimp_nodes

#include <MOM_memory.h>

! Particle status flags
integer, parameter :: ALIVE   = 1 !< Active particle
integer, parameter :: DEAD    = 0 !< Inactive (to be compacted out)
integer, parameter :: LEAVING = 2 !< Particle about to migrate to another PE

! Transfer type flags
integer, parameter :: TRANSFER_FLIP = 0 !< FLIP velocity transfer
integer, parameter :: TRANSFER_APIC = 1 !< APIC (Affine Particle-in-Cell) transfer

! Shape function basis flags
integer, parameter :: BASIS_GIMP  = 1 !< GIMP extended-support shape functions
integer, parameter :: BASIS_SMPM  = 2 !< Standard bilinear (tent) shape functions

integer, parameter :: L_CORNERS   = 1 !< Update corners to get new Lx and Ly
integer, parameter :: L_STRAIN    = 2 !< Calculate Lx and Ly from strain

!> Control structure for the sMPM ice shelf component.
!! All particle arrays are 1-D, indexed 1..n_particles (active) up to n_alloc (allocated).
!! Cell-particle connectivity is maintained in cell_start/cell_count/cell_list.
type, public :: MPM_CS

  ! --- Particle arrays (size n_alloc) ---

  real,    allocatable :: xi(:)       !< Cell-local x-coordinate ∈ [-1,1] [nondim]
  real,    allocatable :: eta(:)      !< Cell-local y-coordinate ∈ [-1,1] [nondim]
  integer, allocatable :: ci(:)       !< T-grid i-index of owning cell
  integer, allocatable :: cj(:)       !< T-grid j-index of owning cell

  real,    allocatable :: up(:)       !< Particle zonal velocity [L T-1 ~> m s-1]
  real,    allocatable :: vp(:)       !< Particle meridional velocity [L T-1 ~> m s-1]
  real,    allocatable :: up_g(:)     !< Previous grid velocity at particle (FLIP) [L T-1 ~> m s-1]
  real,    allocatable :: vp_g(:)     !< Previous grid velocity at particle (FLIP) [L T-1 ~> m s-1]

  real,    allocatable :: H(:)        !< Ice thickness [Z ~> m]
  real,    allocatable :: PVolume(:)  !< Particle area [L2 ~> m2]
  real,    allocatable :: GVolume(:)  !< Original (reference) particle area [L2 ~> m2]
  real,    allocatable :: Lx(:)       !< Particle full length in x [L ~> m]
  real,    allocatable :: Ly(:)       !< Particle full length in y [L ~> m]
  real,    allocatable :: Lx0(:)      !< Original full length in x [L ~> m]
  real,    allocatable :: Ly0(:)      !< Original full length in y [L ~> m]
  real,    allocatable :: strain_x(:) !< Accumulated longitudinal strain ε_xx [nondim]
  real,    allocatable :: strain_y(:) !< Accumulated longitudinal strain ε_yy [nondim]

  !> Deformation gradient F stored as [F11, F12, F21, F22] (row-major 2x2)
  real,    allocatable :: Fdef(:,:)   !< Deformation gradient [nondim]; shape (4, n_alloc)

  real,    allocatable :: GradVel(:,:) !< Velocity gradient [du/dx, dv/dy, du/dy, dv/dx]
                                       !! [T-1 ~> s-1]; shape (4, n_alloc)
  real,    allocatable :: GradZs(:,:) !< Surface-elevation gradient [dZs/dx, dZs/dy]
                                       !! [nondim]; shape (2, n_alloc)
  real,    allocatable :: GradH(:,:)  !< Thickness gradient [dH/dx, dH/dy]
                                       !! [Z L-1 ~> m m-1]; shape (2, n_alloc)

  real,    allocatable :: AGLen(:)    !< Glen's-law rate factor A [Pa-n s-1 ~> R-n L2n T-(1+n)]
                                       !! stored as AGlen_visc equivalent
  real,    allocatable :: eta_visc(:) !< Per-particle depth-integrated effective viscosity
                                       !! eta = 0.5*H*A^(-1/n)*eps_e^((1-n)/n) [R L4 Z T-1 ~> kg m2 s-1]
  real,    allocatable :: newton_visc_factor(:) !< Newton viscosity correction factor
                                       !! = (0.5*(1/n-1)/eps_e2)*eta_visc [R L4 Z T-1 ~> kg m2 s-1]
  real,    allocatable :: basal_trac_p(:)  !< Per-particle linearized basal drag coefficient
                                       !! = beta*H*PVolume [R L2 T-1 ~> kg s-1]; zero for floating
  real,    allocatable :: newton_drag_coef_p(:) !< Per-particle Newton basal drag correction coefficient
                                       !! = (m-1)*basal_trac_p/|u_mid|^2 [R L2 Z-2 T ~> kg m-2 s]; zero for floating
  real,    allocatable :: up_mid(:)    !< Per-particle midpoint u-velocity from Picard iterate [L T-1 ~> m s-1]
  real,    allocatable :: vp_mid(:)    !< Per-particle midpoint v-velocity from Picard iterate [L T-1 ~> m s-1]
  real,    allocatable :: Bp(:,:)     !< APIC affine velocity matrix [L T-1 ~> m s-1]
                                       !! Bp(1,:)=Buu, Bp(2,:)=Buv, Bp(3,:)=Bvu, Bp(4,:)=Bvv
                                       !! Bvv = Σ_I N_I * v_I * (y_I - y_p); shape (4, n_alloc)
  real,    allocatable :: smb(:)      !< Surface mass balance [Z T-1 ~> m s-1]

  integer, allocatable :: status(:)   !< ALIVE / DEAD / LEAVING
  integer(kind=8), allocatable :: pid(:)  !< Globally unique particle ID

  ! --- Cell-particle connectivity (T-grid size, isd:ied × jsd:jed) ---

  integer, allocatable :: cell_count(:,:)  !< Number of active particles in cell (i,j)
  integer, allocatable :: cell_start(:,:)  !< First index in cell_list for cell (i,j)
  integer, allocatable :: cell_list(:)     !< Particle indices in cell order (length n_alloc)

  ! --- B-grid nodal fields (IsdB:IedB × JsdB:JedB) ---

  real, allocatable :: H_node(:,:)    !< Ice thickness interpolated to Bu corners [Z ~> m]
  real, allocatable :: Zs_node(:,:)   !< Surface elevation at Bu corners [Z ~> m]
  real, allocatable :: bed_node(:,:)  !< Bed elevation at Bu corners [Z ~> m]
  real, allocatable :: h_node_mask(:,:)    !< MPM Dirichlet thickness mask at Bu corners [nondim]
                                            !! 1 = prescribed (Dirichlet) node, 0 = free node
  real, allocatable :: h_node_bdry_val(:,:) !< Prescribed ice thickness at Dirichlet Bu corners [Z ~> m]

  ! --- T-grid fields for P2G bookkeeping ---

  real, allocatable :: area_shelf_h(:,:) !< Sum of PVolume_p over particles in cell [L2 ~> m2]
  real, allocatable :: reweight(:,:)     !< areaT / max(area_shelf_h, eps) [nondim]
  logical, allocatable :: is_front_cell(:,:) !< True if cell is ice-covered and adjacent to empty cell

  ! --- Counters and parameters ---

  integer :: n_active  !< Number of active particles on this PE
  integer :: n_alloc   !< Allocated size of particle arrays

  integer :: n_per_cell    !< Number of particles per cell at initialisation (e.g. 4)
  real :: split_factor     !< Split when Lx or Ly > split_factor * Lx0 [nondim]
  real :: flip_alpha       !< FLIP blending coefficient (1=pure FLIP, 0=pure PIC) [nondim]
  integer :: transfer_type !< Particle-grid transfer scheme: 0=FLIP (default), 1=APIC
  integer :: basis_type   !< Shape function basis: BASIS_GIMP (1) or BASIS_SMPM (2)
  integer :: length_update_type   !< Lengths update method: L_CORNERS (1) or L_STRAIN (2)
  real :: H_min         !< Minimum ice thickness at a particle [Z ~> m]
  real :: area_tiny     !< Small area used to avoid division by zero [L2 ~> m2]

  integer(kind=8) :: pid_counter !< Next particle ID to assign on this PE (PE-offset ensures uniqueness)

  logical :: use_constant_AGlen  !< If true, use AGlen_const for all particles
  real    :: AGlen_const         !< Constant Glen A [same units as CS%AGLen]

  real :: n_glen          !< Glen exponent [nondim]
  real :: eps_glen_min    !< Minimum effective strain rate [T-1 ~> s-1]
  real :: density_ice     !< Ice density [R ~> kg m-3]
  real :: density_ocean   !< Ocean density [R ~> kg m-3]
  real :: g_Earth         !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]

  real    :: model_time_total = 0.0 !< Cumulative model time elapsed [T ~> s]; used for VTU filename step number
  real    :: vtu_dt = 0.0          !< VTU output interval [T ~> s]; 0 = off
  real    :: vtu_time_accum = 0.0  !< Time accumulated since last VTU write [T ~> s]
  character(len=200) :: vtu_dir = "VTU"  !< Directory for VTU output files

  ! --- 1D flow-band test (Huth et al. 2021 Section 5.1) ---
  logical :: do_1d_test = .false.  !< If true, enforce analytical H/u in hmask=3 cells
  real    :: test1d_inflow_x_km    !< x-position of inflow boundary [km, geoLonT units]
  real    :: test1d_H0_si          !< Inflow boundary thickness [m] (SI, unscaled)
  real    :: test1d_v0_si          !< Inflow boundary velocity [m s-1] (SI, unscaled)
  real    :: m_to_Z                !< Unit conversion: m -> Z [Z m-1 ~> 1]
  real    :: m_s_to_L_T            !< Unit conversion: m s-1 -> L T-1 [L s T-1 m-1 ~> 1]
  real    :: L_to_m                !< Unit conversion: L -> m [m L-1 ~> 1]

  logical :: do_reweight = .false. !< If true, apply areaT/area_shelf_h reweighting in stiffness assembly
                                    !! (appropriate for sMPM; for GIMP use raw PVolume, i.e. do_reweight=F)

  logical :: initialized = .false. !< True after MPM_init has run

end type MPM_CS

contains

! ============================================================
!> Allocate and initialise the MPM control structure.
!! Seeds particles from ISS%h_shelf (T-grid) and sets up connectivity.
subroutine MPM_init(MPM, ISS, G, US, param_file, n_glen, eps_glen_min, AGLen_visc_ref, &
                    density_ice, density_ocean, g_Earth)
  type(MPM_CS),          intent(inout) :: MPM         !< MPM control structure to initialise
  type(ice_shelf_state), intent(in)    :: ISS         !< Ice shelf state (h_shelf, hmask)
  type(ocean_grid_type), intent(in)    :: G           !< Ocean grid
  type(unit_scale_type), intent(in)    :: US          !< Unit scaling factors
  type(param_file_type), intent(in)    :: param_file  !< Run-time parameter file
  real,                  intent(in)    :: n_glen      !< Glen exponent [nondim]
  real,                  intent(in)    :: eps_glen_min !< Minimum effective strain rate [T-1 ~> s-1]
  real,                  intent(in)    :: AGLen_visc_ref !< Reference Glen A (same units as CS)
  real,                  intent(in)    :: density_ice   !< Ice density [R ~> kg m-3]
  real,                  intent(in)    :: density_ocean !< Ocean density [R ~> kg m-3]
  real,                  intent(in)    :: g_Earth       !< Gravitational acceleration [L2 Z-1 T-2 ~> m s-2]

  ! Local variables
  character(len=40) :: mdl = "MOM_ice_shelf_MPM"
  character(len=9)  :: transfer_str, basis_str, length_update_str
  integer :: i, j, p, pp, ip, jp, n_sq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: pe_rank, npes
  real :: xi0, eta0, dxi, deta, dx, dy, cell_area
  real :: H_p, Vol_p, Lx_p, Ly_p
  real :: dx_km, x_km, x_rel_m, H_si, u_si
  real :: Q0_si, drive_si, alpha_si, m1_si, m2_si
  integer(kind=8) :: id_offset

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  ! --- Read MPM parameters ---
  call get_param(param_file, mdl, "MPM_N_PER_CELL", MPM%n_per_cell, &
                 "Number of material points to seed per ice-covered cell "//&
                 "(must be a perfect square, e.g. 4 or 9).", &
                 default=4, do_not_log=.false.)
  n_sq = nint(sqrt(real(MPM%n_per_cell)))
  if (n_sq*n_sq /= MPM%n_per_cell) &
    call MOM_error(FATAL, "MPM_init: MPM_N_PER_CELL must be a perfect square.")

  call get_param(param_file, mdl, "MPM_SPLIT_FACTOR", MPM%split_factor, &
                 "Split particle when its length exceeds this factor times "//&
                 "original length (Huth et al. 2021 use 1.5).", &
                 units="nondim", default=1.5)
  call get_param(param_file, mdl, "MPM_FLIP_ALPHA", MPM%flip_alpha, &
                 "FLIP blending coefficient; 1 = pure FLIP (no damping), "//&
                 "0 = pure PIC. Ignored when MPM_TRANSFER = APIC.", &
                 units="nondim", default=1.0)

  ! Particle-grid transfer scheme
  call get_param(param_file, mdl, "MPM_TRANSFER", transfer_str, &
                 "Particle-grid velocity transfer scheme. FLIP uses the "//&
                 "FLIP/PIC blended update. APIC (Affine Particle-in-Cell) "//&
                 "uses a local affine correction that eliminates angular "//&
                 "momentum error without artificial viscosity.", &
                 default="APIC")
  if (trim(transfer_str) == "APIC") then
    MPM%transfer_type = TRANSFER_APIC
  else
    MPM%transfer_type = TRANSFER_FLIP
  endif

  call get_param(param_file, mdl, "MPM_REWEIGHT", MPM%do_reweight, &
                 "If true, multiply particle volumes by areaT/area_shelf_h when assembling "//&
                 "the stiffness matrix (sMPM-style reweighting to fill cells). "//&
                 "For GIMP, leave false so raw PVolume is used.", &
                 default=.false.)

  call get_param(param_file, mdl, "MPM_BASIS", basis_str, &
                 "Shape function basis for material points. 'GIMP' (default) uses "//&
                 "generalized interpolation shape functions with extended support "//&
                 "proportional to particle size (Bardenhagen & Kober 2004; Huth "//&
                 "et al. 2021 §4.1). 'sMPM' uses standard bilinear tent functions "//&
                 "with compact nodal support (4 nodes per particle).", &
                 default="GIMP")
  select case (trim(basis_str))
    case ("sMPM")
      MPM%basis_type = BASIS_SMPM
    case ("GIMP")
      MPM%basis_type = BASIS_GIMP
    case default
      MPM%basis_type = BASIS_SMPM
  end select

  call get_param(param_file, mdl, "LENGTH_UPDATE", length_update_str, &
                 "Type of update for material point lengths. 'L_CORNERS' to use the "//&
                 "corner/stretch updating scheme and 'L_STRAIN' to update by tracking "//&
                 "strain (see Huth et al 2021). L_CORNERS can actually by used for all "//&
                 "shape functions, not just GIMP.", default="L_CORNERS")
  select case (trim(length_update_str))
    case ("L_CORNERS")
      MPM%length_update_type = L_CORNERS
    case ("L_STRAIN")
      MPM%length_update_type = L_STRAIN
    case default
      MPM%length_update_type = L_CORNERS
  end select

  call get_param(param_file, mdl, "MPM_H_MIN", MPM%H_min, &
                 "Minimum ice thickness per particle.", &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "MPM_USE_CONSTANT_AGLEN", MPM%use_constant_AGlen, &
                 "If true, use the reference Glen A for all particles.", &
                 default=.true.)
  MPM%AGlen_const = AGLen_visc_ref
  MPM%n_glen        = n_glen
  MPM%eps_glen_min  = eps_glen_min
  MPM%density_ice   = density_ice
  MPM%density_ocean = density_ocean
  MPM%g_Earth       = g_Earth

  call get_param(param_file, mdl, "MPM_VTU_DIR", MPM%vtu_dir, &
                 "Directory for VTU particle visualization files.", &
                 default="VTU")
  call get_param(param_file, mdl, "MPM_VTU_DT", MPM%vtu_dt, &
                 "Time interval for VTU particle output, in days. "//&
                 "0 disables VTU output.", units="days", default=0.0, &
                 scale=86400.0*US%s_to_T)

  ! --- 1D flow-band test parameters (Huth et al. 2021 Section 5.1) ---
  call get_param(param_file, mdl, "MPM_1D_TEST", MPM%do_1d_test, &
                 "If true, enforce analytical Weertman thickness and velocity on "//&
                 "particles in hmask=3 cells for the 1D flow-band verification test.", &
                 default=.false.)
  if (MPM%do_1d_test) then
    call get_param(param_file, mdl, "MPM_1D_INFLOW_X", MPM%test1d_inflow_x_km, &
                   "x-position of the inflow boundary in km (geoLon units). "//&
                   "Extension domain is west of this; active domain is east.", &
                   units="km", default=0.0)
    call get_param(param_file, mdl, "MPM_1D_H0", MPM%test1d_H0_si, &
                   "Inflow boundary ice thickness for 1D test.", &
                   units="m", default=600.0)
    call get_param(param_file, mdl, "MPM_1D_V0", MPM%test1d_v0_si, &
                   "Inflow boundary ice velocity for 1D test.", &
                   units="m yr-1", default=300.0)
    ! Convert v0 from m/yr to m/s (stored unscaled in SI for analytical formula)
    MPM%test1d_v0_si = MPM%test1d_v0_si / 31536000.0 !31556926.0
  endif
  ! Store unit conversion factors for the enforce routine
  MPM%m_to_Z = US%m_to_Z
  MPM%m_s_to_L_T = US%m_s_to_L_T
  MPM%L_to_m = 1.0 / US%m_to_L

  MPM%area_tiny = 1.0e-6 * (US%m_to_L)**2  ! 1 mm^2 in L^2 units

  ! --- Count cells to seed ---
  ! Upper bound: all computational cells could be ice-covered
  MPM%n_alloc = MPM%n_per_cell * (G%iec-G%isc+1) * (G%jec-G%jsc+1) * 2  ! factor 2 headroom for splits
  MPM%n_active = 0

  call MPM_alloc_arrays(MPM, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB)

  ! --- Build PE-unique ID offset ---
  ! Use PE rank × large stride so IDs never collide across PEs
  pe_rank = 0  ! will be overwritten if MPI is active
  npes = 1
  ! Simple approach: pid_counter starts at pe_rank * 10^12
  ! In practice, the ice shelf domain is << 10^12 particles per PE
  MPM%pid_counter = 0  ! IDs assigned sequentially; overlap checked at restart

  ! --- Seed particles from ISS%h_shelf ---
  dxi  = 2.0 / n_sq   ! spacing in local xi  ∈ [-1,1]
  deta = 2.0 / n_sq   ! spacing in local eta ∈ [-1,1]

  MPM%cell_count(:,:) = 0

  ! Seed only the computational domain (isc:iec, jsc:jec); halo cells belong to other PEs
  p = 0
  do j = G%jsc, G%jec
    do i = G%isc, G%iec
      if (ISS%hmask(i,j) /= 1 .and. ISS%hmask(i,j) /= 3) cycle
      H_p = max(ISS%h_shelf(i,j), MPM%H_min)
      dx = G%dxT(i,j)  ! cell width [L]
      dy = G%dyT(i,j)  ! cell height [L]
      cell_area = G%areaT(i,j)  ! [L2]
      Vol_p = cell_area / real(MPM%n_per_cell)  ! particle area [L2]
      Lx_p  = dx / real(n_sq)   ! particle full length in x [L]
      Ly_p  = dy / real(n_sq)   ! particle full length in y [L]

      do jp = 1, n_sq
        do ip = 1, n_sq
          p = p + 1
          if (p > MPM%n_alloc) then
            ! Reallocate (shouldn't happen with initial 2× headroom)
            call MPM_grow_arrays(MPM, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB)
          endif
          ! Place particles on a sub-cell grid, centred in each sub-cell
          xi0  = -1.0 + (real(ip) - 0.5) * dxi
          eta0 = -1.0 + (real(jp) - 0.5) * deta
          MPM%xi(p)  = xi0 ; MPM%eta(p) = eta0
          MPM%ci(p)  = i   ; MPM%cj(p)  = j
          MPM%up(p)  = 0.0 ; MPM%vp(p)  = 0.0
          MPM%up_g(p) = 0.0 ; MPM%vp_g(p) = 0.0
          MPM%Bp(:,p) = 0.0
          MPM%H(p) = H_p
          MPM%PVolume(p) = Vol_p ; MPM%GVolume(p) = Vol_p
          MPM%Lx(p) = Lx_p ; MPM%Ly(p) = Ly_p
          MPM%Lx0(p) = Lx_p ; MPM%Ly0(p) = Ly_p
          MPM%strain_x(p) = 0.0 ; MPM%strain_y(p) = 0.0
          ! Identity deformation gradient
          MPM%Fdef(1,p) = 1.0 ; MPM%Fdef(2,p) = 0.0
          MPM%Fdef(3,p) = 0.0 ; MPM%Fdef(4,p) = 1.0
          MPM%GradVel(1,p) = 0.0 ; MPM%GradVel(2,p) = 0.0
          MPM%GradVel(3,p) = 0.0 ; MPM%GradVel(4,p) = 0.0
          MPM%GradZs(1,p)  = 0.0 ; MPM%GradZs(2,p)  = 0.0
          MPM%GradH(1,p)   = 0.0 ; MPM%GradH(2,p)   = 0.0
          if (MPM%use_constant_AGlen) then
            MPM%AGLen(p) = MPM%AGlen_const
          endif
          MPM%smb(p)    = 0.0
          MPM%status(p) = ALIVE
          MPM%pid(p)    = MPM%pid_counter
          MPM%pid_counter = MPM%pid_counter + 1
          MPM%cell_count(i,j) = MPM%cell_count(i,j) + 1
        enddo
      enddo
    enddo
  enddo

  MPM%n_active = p

  ! --- Apply analytical IC to active-domain particles (1D test only) ---
  ! Active-domain particles (hmask==1) start with H and u matching the
  ! Weertman steady-state profile (same formula as MPM_enforce_1d_test).
  ! This avoids a spin-up transient from zero-velocity initial conditions.
  if (MPM%do_1d_test) then
    Q0_si    = MPM%test1d_H0_si * MPM%test1d_v0_si
    drive_si = 910.0 * 9.81 * (1.0 - 910.0/1028.0) / 4.0
    alpha_si = MPM%AGlen_const * drive_si**MPM%n_glen
    m1_si    = 4.0 * alpha_si / Q0_si
    m2_si    = 1.0 / (MPM%test1d_H0_si**4)
    do p = 1, MPM%n_active
      if (MPM%status(p) /= ALIVE) cycle
      i = MPM%ci(p) ; j = MPM%cj(p)
      if (ISS%hmask(i,j) < 0 ) cycle
      dx_km   = G%dxT(i,j) * MPM%L_to_m * 1.0e-3
      x_km    = G%geoLonT(i,j) + MPM%xi(p) * 0.5 * dx_km
      x_rel_m = (x_km - MPM%test1d_inflow_x_km) * 1.0e3
      if (x_rel_m > 0.0) then
        H_si      = (m1_si * x_rel_m + m2_si)**(-0.25)
        u_si      = Q0_si / H_si
        MPM%H(p)  = H_si * MPM%m_to_Z
        MPM%up(p) = u_si * MPM%m_s_to_L_T
        MPM%vp(p) = 0.0
      else
        MPM%H(p)  = MPM%test1d_H0_si * MPM%m_to_Z
        MPM%up(p) = MPM%test1d_v0_si * MPM%m_s_to_L_T
        MPM%vp(p) = 0.0
      endif
    enddo
  endif

  ! --- Build cell_start / cell_list ---
  call MPM_build_cell_list(MPM, isd, ied, jsd, jed)

  ! --- Initialise B-grid fields ---
  MPM%H_node(:,:) = 0.0 ; MPM%Zs_node(:,:) = 0.0
  MPM%bed_node(:,:) = 0.0
  MPM%h_node_mask(:,:) = 0.0 ; MPM%h_node_bdry_val(:,:) = 0.0
  MPM%area_shelf_h(:,:) = 0.0 ; MPM%reweight(:,:) = 0.0

  MPM%initialized = .true.

  if (is_root_pe()) then
    call MOM_mesg("MPM_init: seeded "//trim(int_to_str(MPM%n_active))//" particles.")
  endif

end subroutine MPM_init

! ============================================================
!> Allocate all MPM arrays.
subroutine MPM_alloc_arrays(MPM, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB)
  type(MPM_CS), intent(inout) :: MPM
  integer, intent(in) :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: n

  n = MPM%n_alloc

  allocate(MPM%xi(n),        source=0.0)
  allocate(MPM%eta(n),       source=0.0)
  allocate(MPM%ci(n),        source=0)
  allocate(MPM%cj(n),        source=0)
  allocate(MPM%up(n),        source=0.0)
  allocate(MPM%vp(n),        source=0.0)
  allocate(MPM%up_g(n),      source=0.0)
  allocate(MPM%vp_g(n),      source=0.0)
  allocate(MPM%H(n),         source=0.0)
  allocate(MPM%PVolume(n),   source=0.0)
  allocate(MPM%GVolume(n),   source=0.0)
  allocate(MPM%Lx(n),        source=0.0)
  allocate(MPM%Ly(n),        source=0.0)
  allocate(MPM%Lx0(n),       source=0.0)
  allocate(MPM%Ly0(n),       source=0.0)
  allocate(MPM%strain_x(n),  source=0.0)
  allocate(MPM%strain_y(n),  source=0.0)
  allocate(MPM%Fdef(4,n),    source=0.0)
  allocate(MPM%GradVel(4,n), source=0.0)
  allocate(MPM%GradZs(2,n),  source=0.0)
  allocate(MPM%GradH(2,n),   source=0.0)
  allocate(MPM%AGLen(n),               source=0.0)
  allocate(MPM%eta_visc(n),            source=0.0)
  allocate(MPM%newton_visc_factor(n),  source=0.0)
  allocate(MPM%basal_trac_p(n),        source=0.0)
  allocate(MPM%newton_drag_coef_p(n),  source=0.0)
  allocate(MPM%up_mid(n),              source=0.0)
  allocate(MPM%vp_mid(n),              source=0.0)
  allocate(MPM%Bp(4,n),                source=0.0)
  allocate(MPM%smb(n),                 source=0.0)
  allocate(MPM%status(n),    source=DEAD)
  allocate(MPM%pid(n),       source=0_8)
  allocate(MPM%cell_list(n), source=0)

  allocate(MPM%cell_count(isd:ied, jsd:jed), source=0)
  allocate(MPM%cell_start(isd:ied, jsd:jed), source=0)

  allocate(MPM%H_node(IsdB:IedB, JsdB:JedB),        source=0.0)
  allocate(MPM%Zs_node(IsdB:IedB, JsdB:JedB),       source=0.0)
  allocate(MPM%bed_node(IsdB:IedB, JsdB:JedB),       source=0.0)
  allocate(MPM%h_node_mask(IsdB:IedB, JsdB:JedB),    source=0.0)
  allocate(MPM%h_node_bdry_val(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(MPM%area_shelf_h(isd:ied, jsd:jed),  source=0.0)
  allocate(MPM%reweight(isd:ied, jsd:jed),       source=0.0)
  allocate(MPM%is_front_cell(isd:ied, jsd:jed))
  MPM%is_front_cell(:,:) = .false.

end subroutine MPM_alloc_arrays

! ============================================================
!> Grow particle arrays by 50% when near capacity.
subroutine MPM_grow_arrays(MPM, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB)
  type(MPM_CS), intent(inout) :: MPM
  integer, intent(in) :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB

  integer :: n_new, n_old
  real, allocatable :: tmp_r1(:), tmp_r2(:,:)
  integer, allocatable :: tmp_i(:)
  integer(kind=8), allocatable :: tmp_i8(:)

  n_old = MPM%n_alloc
  n_new = n_old + n_old/2 + 1
  MPM%n_alloc = n_new

  call grow_real_1d(MPM%xi, n_new) ; call grow_real_1d(MPM%eta, n_new)
  call grow_int_1d(MPM%ci, n_new)  ; call grow_int_1d(MPM%cj, n_new)
  call grow_real_1d(MPM%up, n_new)  ; call grow_real_1d(MPM%vp, n_new)
  call grow_real_1d(MPM%up_g, n_new) ; call grow_real_1d(MPM%vp_g, n_new)
  call grow_real_1d(MPM%H, n_new)
  call grow_real_1d(MPM%PVolume, n_new) ; call grow_real_1d(MPM%GVolume, n_new)
  call grow_real_1d(MPM%Lx, n_new)  ; call grow_real_1d(MPM%Ly, n_new)
  call grow_real_1d(MPM%Lx0, n_new) ; call grow_real_1d(MPM%Ly0, n_new)
  call grow_real_1d(MPM%strain_x, n_new) ; call grow_real_1d(MPM%strain_y, n_new)
  call grow_real_2d(MPM%Fdef, 4, n_new)
  call grow_real_2d(MPM%GradVel, 4, n_new)
  call grow_real_2d(MPM%GradZs, 2, n_new)
  call grow_real_2d(MPM%GradH,  2, n_new)
  call grow_real_1d(MPM%AGLen, n_new)
  call grow_real_1d(MPM%eta_visc, n_new)
  call grow_real_1d(MPM%newton_visc_factor, n_new)
  call grow_real_1d(MPM%basal_trac_p, n_new)
  call grow_real_1d(MPM%newton_drag_coef_p, n_new)
  call grow_real_1d(MPM%up_mid, n_new)
  call grow_real_1d(MPM%vp_mid, n_new)
  call grow_real_2d(MPM%Bp, 4, n_new)
  call grow_real_1d(MPM%smb, n_new)
  call grow_int_1d(MPM%status, n_new)
  call grow_int8_1d(MPM%pid, n_new)
  call grow_int_1d(MPM%cell_list, n_new)

end subroutine MPM_grow_arrays

! ============================================================
!> Rebuild the cell_start and cell_list arrays from cell_count.
subroutine MPM_build_cell_list(MPM, isd, ied, jsd, jed)
  type(MPM_CS), intent(inout) :: MPM
  integer, intent(in) :: isd, ied, jsd, jed

  integer :: i, j, p, offset
  integer, allocatable :: fill(:,:)

  allocate(fill(isd:ied, jsd:jed), source=0)

  ! Compute cell_start by prefix sum over cell_count
  offset = 0
  do j = jsd, jed
    do i = isd, ied
      MPM%cell_start(i,j) = offset + 1
      offset = offset + MPM%cell_count(i,j)
    enddo
  enddo

  ! Place particle indices into cell_list
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    i = MPM%ci(p) ; j = MPM%cj(p)
    fill(i,j) = fill(i,j) + 1
    MPM%cell_list(MPM%cell_start(i,j) + fill(i,j) - 1) = p
  enddo

  deallocate(fill)

end subroutine MPM_build_cell_list

! ============================================================
!> Particles → Grid (P2G): update ISS%h_shelf, ISS%hmask, MPM%H_node, MPM%Zs_node.
!!
!! B-grid H_node is computed by GIMP-weighted scatter from particles.
!! h_shelf for each T-grid cell is then the average of its 4 B-grid corner H_nodes,
!! which activates any cell overlapped by a GIMP particle domain (Huth et al. 2021):
!! "any element that a GIMP domain overlaps should be an active element."
!! The T-grid reweight factor (areaT / Σ PVolume) is retained for stiffness quadrature.
subroutine MPM_P2G(MPM, ISS, G)
  type(MPM_CS),          intent(inout) :: MPM  !< MPM control structure
  type(ice_shelf_state), intent(inout) :: ISS  !< Ice shelf state (h_shelf, hmask updated here)
  type(ocean_grid_type), intent(in)    :: G    !< Ocean grid

  integer :: i, j, p, ic, jc, I_Bu, J_Bu, k_n, n_gimp
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real, allocatable :: H_vol(:,:)  ! Σ H(p)*PVolume(p) per cell [Z L2]
  real, allocatable :: H_wgt(:,:)  ! Σ PVolume(p) per cell [L2]
  real, allocatable :: Hn_vol(:,:) ! Σ N(k)*H(p)*w_eff per Bu corner [Z L2]
  real, allocatable :: Hn_wgt(:,:) ! Σ N(k)*w_eff per Bu corner [L2]
  real :: N_g(9), dNdx_g(9), dNdy_g(9)
  integer :: di_g(9), dj_g(9)
  real :: w_eff, dens_ratio, dx, dy, H_cell

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  allocate(H_vol(isd:ied, jsd:jed), source=0.0)
  allocate(H_wgt(isd:ied, jsd:jed), source=0.0)
  allocate(Hn_vol(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(Hn_wgt(IsdB:IedB, JsdB:JedB), source=0.0)

  do j = jsd, jed ; do i = isd, ied
    if (ISS%hmask(i,j) == 1) then
      ISS%hmask(i,j) = 0
      ISS%h_shelf(i,j) = 0.0
    endif
  enddo; enddo

  ! --- Accumulate per-cell particle H and area ---
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    H_wgt(ic,jc) = H_wgt(ic,jc) + MPM%PVolume(p)
    H_vol(ic,jc) = H_vol(ic,jc) + MPM%H(p) * MPM%PVolume(p)
    if (ISS%hmask(ic,jc) == 0) ISS%hmask(ic,jc) = 1   ! cell activated by GIMP domain overlap
  enddo

  ! --- Update area_shelf_h and reweight ---
  MPM%area_shelf_h(:,:) = H_wgt(:,:)
  if (MPM%do_reweight) then
    ! sMPM-style: scale particle volumes so they exactly fill each cell
    do j = jsd, jed ; do i = isd, ied
      if (MPM%area_shelf_h(i,j) > MPM%area_tiny) then
        MPM%reweight(i,j) = G%areaT(i,j) / MPM%area_shelf_h(i,j)
      else
        MPM%reweight(i,j) = 0.0
      endif
    enddo ; enddo
  else
    ! GIMP: use raw PVolume; reweight = 1 everywhere
    MPM%reweight(:,:) = 1.0
  endif

  ! --- B-grid P2G: interpolate particle H to Bu corners using GIMP shape functions ---
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    dx = G%dxT(ic,jc) ; dy = G%dyT(ic,jc)
    if (MPM%basis_type == BASIS_GIMP) then
      call gimp_nodes(MPM%xi(p), MPM%eta(p), 0.5*MPM%Lx(p), 0.5*MPM%Ly(p), dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    else
      call smpm_nodes(MPM%xi(p), MPM%eta(p), dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    endif
    ! Eq 29 (Huth 2021): use raw PVolume, not reweighted area
    ! reweight is only used in CG_action_MPM stiffness assembly
    w_eff = MPM%PVolume(p)

    ! Scatter to all contributing Bu corners (standard + GIMP extended)
    do k_n = 1, n_gimp
      I_Bu = ic - 1 + di_g(k_n)
      J_Bu = jc - 1 + dj_g(k_n)
      if (I_Bu < IsdB .or. I_Bu > IedB) cycle
      if (J_Bu < JsdB .or. J_Bu > JedB) cycle
      Hn_wgt(I_Bu, J_Bu) = Hn_wgt(I_Bu, J_Bu) + N_g(k_n) * w_eff
      Hn_vol(I_Bu, J_Bu) = Hn_vol(I_Bu, J_Bu) + N_g(k_n) * MPM%H(p) * w_eff
    enddo
  enddo

  ! Halo exchange before dividing
  call pass_var(Hn_wgt, G%domain, position=CORNER)
  call pass_var(Hn_vol, G%domain, position=CORNER)

  ! Compute H_node and Zs_node at B-grid corners. TODO: Add grounded Zs
  dens_ratio = MPM%density_ice / MPM%density_ocean
  do J_Bu = JsdB, JedB ; do I_Bu = IsdB, IedB
    if (Hn_wgt(I_Bu, J_Bu) > MPM%area_tiny) then
      MPM%H_node(I_Bu, J_Bu)  = Hn_vol(I_Bu, J_Bu) / Hn_wgt(I_Bu, J_Bu)
      MPM%Zs_node(I_Bu, J_Bu) = MPM%H_node(I_Bu, J_Bu) * (1.0 - dens_ratio)
    else
      MPM%H_node(I_Bu, J_Bu)  = 0.0
      MPM%Zs_node(I_Bu, J_Bu) = 0.0
    endif
  enddo ; enddo

  ! Enforce prescribed thickness at MPM Dirichlet nodes (h_node_mask==1).
  ! This overrides any particle-based H_node values in the inflow BC region.
  do J_Bu = JsdB, JedB ; do I_Bu = IsdB, IedB
    if (MPM%h_node_mask(I_Bu, J_Bu) == 1.0) then
      MPM%H_node(I_Bu, J_Bu)  = MPM%h_node_bdry_val(I_Bu, J_Bu)
      MPM%Zs_node(I_Bu, J_Bu) = MPM%H_node(I_Bu, J_Bu) * (1.0 - dens_ratio)
    endif
  enddo ; enddo

  ! --- Update ISS%h_shelf and ISS%hmask via H_node corner averages --- THIS IS DONE SO WRONG
  ! GIMP activation criterion (Huth et al. 2021): any element whose domain a
  ! GIMP particle overlaps must be an active element.  H_node was computed with
  ! full GIMP support above, so extended corners adjacent to the ice front are
  ! non-zero even before a particle center crosses the cell boundary.  Averaging
  ! the 4 B-grid corners gives h_shelf for each T-grid cell and naturally
  ! activates cells that lie within a GIMP particle domain.
  ! Only update the computational domain (isc..iec, jsc..jec).  Halo cells must
  ! NOT be activated here: a halo cell adjacent to the calving front has two
  ! valid corners (H_node > 0) and two OOB corners (H_node = 0), giving H_cell > 0
  ! even though no ice exists there.  Activating it would make the last true ice
  ! cell appear to have an active eastern neighbour, suppressing front detection.
  ! After pass_var below, halos are filled with correct values from neighbouring PEs.
  do j = G%jsc, G%jec ; do i = G%isc, G%iec
    H_cell = 0.25 * (MPM%H_node(i-1,j-1) + MPM%H_node(i,j-1) + &
                     MPM%H_node(i-1,j  ) + MPM%H_node(i,j  ))
    if (ISS%hmask(i,j) == 1) then
      ISS%h_shelf(i,j) = max(H_cell, MPM%H_min)
      !if (ISS%hmask(i,j) == 0) ISS%hmask(i,j) = 1   ! cell activated by GIMP domain overlap
    ! else
    !   if (ISS%hmask(i,j) == 1) then
    !     ISS%hmask(i,j) = 0
    !     ISS%h_shelf(i,j) = 0.0
    !   endif
    endif
  enddo ; enddo

  ! Propagate h_shelf and hmask into halo so the SSA solver sees them
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%hmask,   G%domain, complete=.true.)

  deallocate(H_vol, H_wgt, Hn_vol, Hn_wgt)

end subroutine MPM_P2G

! ============================================================
!> Update masks (umask, vmask, hmask front cells) based on current hmask.
!! Called by the dynamics module after MPM_P2G.
subroutine MPM_update_masks(MPM, ISS, G)
  type(MPM_CS),          intent(in)    :: MPM
  type(ice_shelf_state), intent(in)    :: ISS
  type(ocean_grid_type), intent(in)    :: G
  ! Mask updates are handled by update_velocity_masks in MOM_ice_shelf_dynamics;
  ! this routine is a hook for future extensions.
end subroutine MPM_update_masks

! ============================================================
!> Grid → Particles (G2P): interpolate u_shelf and Zs to particles,
!! computing particle velocity (for FLIP) and velocity/Zs gradients.
subroutine MPM_G2P(MPM, G, u_shelf, v_shelf, hmask)
  type(MPM_CS),          intent(inout) :: MPM     !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G       !< Ocean grid
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(in)    :: u_shelf !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(in)    :: v_shelf !< Meridional velocity [L T-1 ~> m s-1]
  real, dimension(G%isd:G%ied, G%jsd:G%jed), &
                         intent(in)    :: hmask   !< Ice-shelf thickness mask (T-grid);
                                                  !! hmask==3 marks prescribed extension cells.

  integer :: p, ic, jc, k_n, n_gimp, I_Bu, J_Bu, IsdB, IedB, JsdB, JedB
  real :: xi_p, eta_p, dx, dy, dx_half, dy_half
  real :: N_g(9), dNdx_g(9), dNdy_g(9)
  integer :: di_g(9), dj_g(9)
  real :: u_I, v_I, Zs_I       ! Velocity and surface elevation at grid node [L T-1], [Z]
  real :: u_new, v_new          ! Interpolated grid velocity at particle [L T-1]
  real :: dx_I_k, dy_I_k        ! Node-to-particle offsets for APIC [L]

  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    xi_p  = MPM%xi(p)
    eta_p = MPM%eta(p)
    dx = G%dxT(ic,jc)
    dy = G%dyT(ic,jc)
    dx_half = 0.5 * dx
    dy_half = 0.5 * dy

    if (MPM%basis_type == BASIS_GIMP) then
      call gimp_nodes(xi_p, eta_p, 0.5*MPM%Lx(p), 0.5*MPM%Ly(p), dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    else
      call smpm_nodes(xi_p, eta_p, dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    endif

    ! G2P gather: accumulate interpolated values from all contributing Bu nodes
    u_new = 0.0 ; v_new = 0.0
    MPM%GradVel(1,p) = 0.0 ; MPM%GradVel(2,p) = 0.0
    MPM%GradVel(3,p) = 0.0 ; MPM%GradVel(4,p) = 0.0
    MPM%GradZs(1,p) = 0.0 ; MPM%GradZs(2,p) = 0.0
    MPM%GradH(1,p)  = 0.0 ; MPM%GradH(2,p)  = 0.0
    if (MPM%transfer_type == TRANSFER_APIC) then
      MPM%Bp(1,p) = 0.0 ; MPM%Bp(2,p) = 0.0
      MPM%Bp(3,p) = 0.0 ; MPM%Bp(4,p) = 0.0
    endif

    do k_n = 1, n_gimp
      I_Bu = ic - 1 + di_g(k_n)
      J_Bu = jc - 1 + dj_g(k_n)
      if (I_Bu < IsdB .or. I_Bu > IedB) cycle
      if (J_Bu < JsdB .or. J_Bu > JedB) cycle

      u_I  = u_shelf(I_Bu, J_Bu)
      v_I  = v_shelf(I_Bu, J_Bu)
      Zs_I = MPM%Zs_node(I_Bu, J_Bu)

      u_new = u_new + N_g(k_n) * u_I
      v_new = v_new + N_g(k_n) * v_I

      ! Velocity gradients ∇u for Glen's law viscosity
      MPM%GradVel(1,p) = MPM%GradVel(1,p) + dNdx_g(k_n) * u_I  ! du/dx
      MPM%GradVel(2,p) = MPM%GradVel(2,p) + dNdy_g(k_n) * v_I  ! dv/dy
      MPM%GradVel(3,p) = MPM%GradVel(3,p) + dNdy_g(k_n) * u_I  ! du/dy
      MPM%GradVel(4,p) = MPM%GradVel(4,p) + dNdx_g(k_n) * v_I  ! dv/dx

      ! Surface elevation gradient for driving stress
      MPM%GradZs(1,p) = MPM%GradZs(1,p) + dNdx_g(k_n) * Zs_I
      MPM%GradZs(2,p) = MPM%GradZs(2,p) + dNdy_g(k_n) * Zs_I

      ! Thickness gradient for child-thickness correction at splitting (Huth 2021 §4.2 eq 44)
      MPM%GradH(1,p) = MPM%GradH(1,p) + dNdx_g(k_n) * MPM%H_node(I_Bu, J_Bu)
      MPM%GradH(2,p) = MPM%GradH(2,p) + dNdy_g(k_n) * MPM%H_node(I_Bu, J_Bu)

      ! APIC: B_p = Σ_I N_I * v_I * (x_I - x_p)^T
      ! Node position in local coords: xi_I = 2*di_g - 1, eta_I = 2*dj_g - 1
      if (MPM%transfer_type == TRANSFER_APIC) then
        dx_I_k = (-1.0 + 2.0*real(di_g(k_n)) - xi_p)  * dx_half
        dy_I_k = (-1.0 + 2.0*real(dj_g(k_n)) - eta_p) * dy_half
        MPM%Bp(1,p) = MPM%Bp(1,p) + N_g(k_n) * u_I * dx_I_k  ! Buu
        MPM%Bp(2,p) = MPM%Bp(2,p) + N_g(k_n) * u_I * dy_I_k  ! Buv
        MPM%Bp(3,p) = MPM%Bp(3,p) + N_g(k_n) * v_I * dx_I_k  ! Bvu
        MPM%Bp(4,p) = MPM%Bp(4,p) + N_g(k_n) * v_I * dy_I_k  ! Bvv
      endif
    enddo

    ! Extension-domain particles (hmask=3) have prescribed velocity set by
    ! MPM_enforce_1d_test. Skip the FLIP/APIC transfer so that the P2G-set
    ! pre-solve grid state cannot corrupt up(p) before enforce runs.
    !if (hmask(ic, jc) == 3.0) cycle

    if (MPM%transfer_type == TRANSFER_APIC) then
      MPM%up(p) = u_new
      MPM%vp(p) = v_new
      ! up_g not used in APIC but update for consistency
      ! MPM%up_g(p) = u_new
      ! MPM%vp_g(p) = v_new
    else
      ! FLIP velocity update (Eq 34, Huth 2021):
      !   v_p^{m+1} = v_p^m + Σ_I (v_I^{m+1} - v_I^m) * N_I
      ! where v_I^m = up_g(p) was set by MPM_store_pre_solve_vel (called before SSA)
      MPM%up(p) = MPM%flip_alpha * (MPM%up(p) + (u_new - MPM%up_g(p))) &
                + (1.0 - MPM%flip_alpha) * u_new
      MPM%vp(p) = MPM%flip_alpha * (MPM%vp(p) + (v_new - MPM%vp_g(p))) &
                + (1.0 - MPM%flip_alpha) * v_new
      ! Save current grid velocity for next FLIP step
      ! MPM%up_g(p) = u_new
      ! MPM%vp_g(p) = v_new
    endif
  enddo

end subroutine MPM_G2P

! ============================================================
!> Lagrangian particle updates: deformation gradient F, thickness H,
!! particle dimensions, and position.  Called once per timestep after G2P.
subroutine MPM_Lagrangian_update(MPM, G, dt, hmask)
  type(MPM_CS),          intent(inout) :: MPM   !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G     !< Ocean grid
  real,                  intent(in)    :: dt    !< Timestep [T ~> s]
  real, dimension(G%isd:G%ied, G%jsd:G%jed), &
                         intent(in)    :: hmask !< Ice-shelf thickness mask; hmask==0 is land

  integer :: p, ic, jc
  real :: ux, uy, vx, vy   ! Velocity gradient components [T-1]
  real :: divv              ! Velocity divergence [T-1]
  real :: F11, F12, F21, F22   ! Deformation gradient (old)
  real :: L11, L12, L21, L22  ! Velocity gradient increment L = I + dt*gradV
  real :: detF              ! Det of updated deformation gradient
  real :: H_new, xi_new, eta_new
  real :: dx, dy
  ! GIMP domain update
  logical :: use_corners
  real :: Lx_half, Ly_half
  real :: mx(4), my(4), tl1, tl2, sterm  ! corner-tracking work arrays
  real :: G11, G12, G22                  ! right Cauchy-Green tensor components
  real :: tr_G, disc_G, lam1, lam2, U1, U2  ! eigenvalues / principal stretches
  real :: q2x, q2y, q1x, q1y, qnorm     ! eigenvectors of G

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle

    ic = MPM%ci(p) ; jc = MPM%cj(p)
    dx = G%dxT(ic, jc)
    dy = G%dyT(ic, jc)

    ! --- Update deformation gradient F^{n+1} = (I + dt*L) · F^n ---
    ux = MPM%GradVel(1,p)  ! ∂u/∂x
    vy = MPM%GradVel(2,p)  ! ∂v/∂y
    uy = MPM%GradVel(3,p)  ! ∂u/∂y
    vx = MPM%GradVel(4,p)  ! ∂v/∂x
    divv = ux + vy

    L11 = 1.0 + dt*ux ; L12 = dt*uy
    L21 = dt*vx       ; L22 = 1.0 + dt*vy

    F11 = MPM%Fdef(1,p) ; F12 = MPM%Fdef(2,p)
    F21 = MPM%Fdef(3,p) ; F22 = MPM%Fdef(4,p)

    MPM%Fdef(1,p) = L11*F11 + L12*F21
    MPM%Fdef(2,p) = L11*F12 + L12*F22
    MPM%Fdef(3,p) = L21*F11 + L22*F21
    MPM%Fdef(4,p) = L21*F12 + L22*F22

    detF = MPM%Fdef(1,p)*MPM%Fdef(4,p) - MPM%Fdef(2,p)*MPM%Fdef(3,p)


    ! --- Update particle domain size and volume ---
    if (MPM%length_update_type == L_STRAIN) then
      MPM%strain_x(p) = MPM%strain_x(p) + dt*ux
      MPM%strain_y(p) = MPM%strain_y(p) + dt*vy

      MPM%PVolume(p) = detF * MPM%GVolume(p)
      MPM%Lx(p) = MPM%Lx0(p) * (1.0 + MPM%strain_x(p))
      MPM%Ly(p) = MPM%Ly0(p) * (1.0 + MPM%strain_y(p))

    elseif (MPM%length_update_type == L_CORNERS) then
      ! Update GIMP support-domain full lengths (Huth et al. 2021 §3.4).
      ! Interior particles: corner-tracking scheme.  Boundary particles: stretch-tensor.
      use_corners = (ic > G%isd .and. ic < G%ied .and. &
                     jc > G%jsd .and. jc < G%jed)
      if (use_corners) use_corners = &
          hmask(ic-1,jc) > 0.5 .and. hmask(ic+1,jc) > 0.5 .and. &
          hmask(ic,jc-1) > 0.5 .and. hmask(ic,jc+1) > 0.5

      if (use_corners) then
        ! Corner-tracking: advect the 4 side-midpoints at (±Lx/2,0) and (0,±Ly/2)
        ! using the local velocity gradient.  Particle advection (up,vp) cancels in
        ! the span, so only gradient terms appear.
        Lx_half = 0.5 * MPM%Lx(p) ; Ly_half = 0.5 * MPM%Ly(p)

        ! x-coords of the 4 midpoints after dt (relative to particle centre):
        mx(1) = -Lx_half * (1.0 + dt*ux)   ! left  midpoint (−Lx/2, 0)
        mx(2) =  dt * uy * Ly_half          ! top   midpoint (0, +Ly/2)
        mx(3) =  Lx_half * (1.0 + dt*ux)   ! right midpoint (+Lx/2, 0)
        mx(4) = -dt * uy * Ly_half          ! bot   midpoint (0, −Ly/2)

        ! y-coords:
        my(1) = -dt * vx * Lx_half          ! left
        my(2) =  Ly_half * (1.0 + dt*vy)   ! top
        my(3) =  dt * vx * Lx_half          ! right
        my(4) = -Ly_half * (1.0 + dt*vy)   ! bottom

        tl1 = maxval(mx) - minval(mx)   ! new full x-length
        tl2 = maxval(my) - minval(my)   ! new full y-length

        ! Volume-correction factor: ensures Lx_new*Ly_new = detF*Lx0*Ly0 = PVolume_new
        if (tl1 > 0.0 .and. tl2 > 0.0) then
          sterm = sqrt(detF * MPM%Lx0(p) * MPM%Ly0(p) / (tl1 * tl2))
          MPM%Lx(p) = tl1 * sterm
          MPM%Ly(p) = tl2 * sterm
        endif

      else
        ! Stretch-tensor fallback (Huth et al. 2021 §3.4 / Bardenhagen & Kober 2004).
        ! Compute right Cauchy-Green tensor G = F^T F (symmetric).
        G11 = MPM%Fdef(1,p)**2 + MPM%Fdef(3,p)**2
        G12 = MPM%Fdef(1,p)*MPM%Fdef(2,p) + MPM%Fdef(3,p)*MPM%Fdef(4,p)
        G22 = MPM%Fdef(2,p)**2 + MPM%Fdef(4,p)**2

        ! 2×2 symmetric eigendecomposition: eigenvalues = principal stretches².
        tr_G   = G11 + G22
        disc_G = sqrt(max(0.25*(G11-G22)**2 + G12**2, 0.0))
        lam1   = 0.5*tr_G - disc_G   ! smaller eigenvalue
        lam2   = 0.5*tr_G + disc_G   ! larger eigenvalue
        U1 = sqrt(max(lam1, 0.0))    ! principal stretches
        U2 = sqrt(max(lam2, 0.0))

        ! Eigenvector for lam2 (larger), normalised.
        if (abs(G12) > 1.0e-30 * abs(tr_G)) then
          q2x = lam2 - G22 ; q2y = G12
        else
          q2x = 1.0 ; q2y = 0.0   ! already diagonal
        endif
        qnorm = sqrt(q2x**2 + q2y**2)
        q2x = q2x/qnorm ; q2y = q2y/qnorm
        q1x = -q2y ; q1y = q2x   ! orthogonal eigenvector for lam1

        ! New full lengths (diagonal of right stretch tensor U in reference basis).
        MPM%Lx(p) = MPM%Lx0(p) * (U1*q1x**2 + U2*q2x**2)
        MPM%Ly(p) = MPM%Ly0(p) * (U1*q1y**2 + U2*q2y**2)
      endif

        MPM%PVolume(p) = MPM%Lx(p) * MPM%Ly(p)
    endif

    ! --- Update ice thickness: dH/dt = -H * div(v) + SMB ---
    H_new = MPM%H(p) * (1.0 - dt * divv) + MPM%smb(p) * dt
    if (H_new < MPM%H_min) H_new = MPM%H_min  ! enforce minimum thickness
    MPM%H(p) = H_new

    ! --- Advect particle position ---
    ! Update local coordinates; migration handles cell crossing
    xi_new  = MPM%xi(p)  + dt * MPM%up(p) / (0.5 * dx)
    eta_new = MPM%eta(p) + dt * MPM%vp(p) / (0.5 * dy)
    MPM%xi(p)  = xi_new
    MPM%eta(p) = eta_new
  enddo

end subroutine MPM_Lagrangian_update

! ============================================================
!> Split particles that have grown too large (Lx or Ly > split_factor * L0).
!! A 2-way split is performed along the longest dimension; a 4-way split when
!! both dimensions exceed the threshold.
subroutine MPM_split(MPM, G)
  type(MPM_CS),          intent(inout) :: MPM  !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G    !< Ocean grid

  integer :: p, n_new, n_before
  integer :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  logical :: split_x, split_y
  real :: half_xi, half_eta, Lx_new, Ly_new, dx, dy
  real :: dxi_child, deta_child

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  n_before = MPM%n_active

  do p = 1, n_before   ! Only loop original particles, not newly created ones
    if (MPM%status(p) /= ALIVE) cycle


    dx = G%dxT(MPM%ci(p), MPM%cj(p))
    dy = G%dyT(MPM%ci(p), MPM%cj(p))

    split_x = (MPM%Lx(p) > MPM%split_factor * (dx/sqrt(real(MPM%n_per_cell))))
    split_y = (MPM%Ly(p) > MPM%split_factor * (dy/sqrt(real(MPM%n_per_cell))))

    if (.not.(split_x .or. split_y)) cycle

    ! Offset in local coordinates to child centre: ±Lx/4 physical = ±(Lx/4)/(dx/2) local
    dxi_child  = MPM%Lx(p) / dx   ! [nondim]; 0.5*dxi_child = Lx/(2*dx) local = Lx/4 physical
    deta_child = MPM%Ly(p) / dy   ! [nondim]

    if (split_x .and. split_y) then
      ! 4-way split: children at ±Lx/4, ±Ly/4 from parent (Huth et al. 2021 §4.2)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 +0.5*dxi_child, +0.5*deta_child, 0.25, dx, dy)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 -0.5*dxi_child, +0.5*deta_child, 0.25, dx, dy)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 +0.5*dxi_child, -0.5*deta_child, 0.25, dx, dy)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 -0.5*dxi_child, -0.5*deta_child, 0.25, dx, dy)
      MPM%status(p) = DEAD  ! parent is replaced by 4 children
    elseif (split_x) then
      ! 2-way split in x: children at ±Lx/2 (Huth et al. 2021 §4.2)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 +0.5*dxi_child, 0.0, 0.5, dx, dy)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 -0.5*dxi_child, 0.0, 0.5, dx, dy)
      MPM%status(p) = DEAD
    else  ! split_y only
      ! 2-way split in y: children at ±Ly/2 (Huth et al. 2021 §4.2)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 0.0, +0.5*deta_child, 0.5, dx, dy)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 0.0, -0.5*deta_child, 0.5, dx, dy)
      MPM%status(p) = DEAD
    endif
  enddo

end subroutine MPM_split

! ============================================================
!> Create a child particle from parent p, offset in local coords.
subroutine split_particle_child(MPM, p_parent, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                  dxi, deta, vol_frac, dx, dy)
  type(MPM_CS), intent(inout) :: MPM
  integer,      intent(in)    :: p_parent
  integer,      intent(in)    :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real,         intent(in)    :: dxi, deta    !< Offset in cell-local coords
  real,         intent(in)    :: vol_frac     !< Fraction of parent volume
  real,         intent(in)    :: dx, dy       !< Physical cell dimensions [L ~> m]

  integer :: p_child

  MPM%n_active = MPM%n_active + 1
  if (MPM%n_active > MPM%n_alloc) then
    call MPM_grow_arrays(MPM, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB)
  endif
  p_child = MPM%n_active

  MPM%xi(p_child)  = MPM%xi(p_parent)  + dxi
  MPM%eta(p_child) = MPM%eta(p_parent) + deta
  MPM%ci(p_child)  = MPM%ci(p_parent)
  MPM%cj(p_child)  = MPM%cj(p_parent)
  MPM%up(p_child)  = MPM%up(p_parent)
  MPM%vp(p_child)  = MPM%vp(p_parent)
  MPM%up_g(p_child) = MPM%up_g(p_parent)
  MPM%vp_g(p_child) = MPM%vp_g(p_parent)

  MPM%PVolume(p_child) = vol_frac * MPM%PVolume(p_parent)
  ! Child full lengths: halved in the split direction(s)
  ! For 4-way: both halved; for 2-way x: Lx halved; for 2-way y: Ly halved
  ! The vol_frac encodes this: 0.25 = 4-way, 0.5 = 2-way
  if (vol_frac < 0.3) then
    ! 4-way split: both dimensions halved
    MPM%Lx(p_child) = 0.5 * MPM%Lx(p_parent)
    MPM%Ly(p_child) = 0.5 * MPM%Ly(p_parent)
    if (MPM%length_update_type == L_CORNERS) then
      MPM%Lx0(p_child) = 0.5 * MPM%Lx0(p_parent)
      MPM%Ly0(p_child) = 0.5 * MPM%Ly0(p_parent)
    else
      MPM%Lx0(p_child) = sqrt(0.5)*MPM%Lx0(p_parent)
      MPM%Ly0(p_child) = MPM%Lx0(p_child)
      MPM%strain_x(p_child) = (MPM%strain_x(p_parent) + 1.0)/ (2.0*sqrt(0.25)) - 1.0
      MPM%strain_y(p_child) = (MPM%strain_y(p_parent) + 1.0)/ (2.0*sqrt(0.25)) - 1.0
    endif
    MPM%H(p_child) = MPM%H(p_parent) + &
      MPM%GradH(1, p_parent) * dxi  * (0.5 * dx) + &
      MPM%GradH(2, p_parent) * deta * (0.5 * dy)
  elseif (abs(deta) < 1.0e-10) then
    ! 2-way x-split: only Lx halved
    MPM%Lx(p_child) = 0.5 * MPM%Lx(p_parent)
    MPM%Ly(p_child) = MPM%Ly(p_parent)
    if (MPM%length_update_type == L_CORNERS) then
      MPM%Lx0(p_child) = 0.5 * MPM%Lx0(p_parent)
      MPM%Ly0(p_child) = MPM%Ly0(p_parent)
    else
      MPM%Lx0(p_child) = sqrt(0.5)*MPM%Lx0(p_parent)
      MPM%Ly0(p_child) = MPM%Lx0(p_child)
      MPM%strain_x(p_child) = (MPM%strain_x(p_parent) + 1.0)/ (2.0*sqrt(0.5)) - 1.0
      MPM%strain_y(p_child) = (MPM%strain_y(p_parent) + 1.0)/ sqrt(0.5) - 1.0
    endif
    MPM%H(p_child) = MPM%H(p_parent) + MPM%GradH(1, p_parent) * dxi * (0.5 * dx)
  else
    ! 2-way y-split: only Ly halved
    MPM%Lx(p_child) = MPM%Lx(p_parent)
    MPM%Ly(p_child) = 0.5 * MPM%Ly(p_parent)
    if (MPM%length_update_type == L_CORNERS) then
      MPM%Lx0(p_child) = MPM%Lx0(p_parent)
      MPM%Ly0(p_child) = 0.5 * MPM%Ly0(p_parent)
    else
      MPM%Lx0(p_child) = sqrt(0.5)*MPM%Lx0(p_parent)
      MPM%Ly0(p_child) = MPM%Lx0(p_child)
      MPM%strain_x(p_child) = (MPM%strain_x(p_parent) + 1.0)/ sqrt(0.5) - 1.0
      MPM%strain_y(p_child) = (MPM%strain_y(p_parent) + 1.0)/ (2.0*sqrt(0.5)) - 1.0
    endif
    MPM%H(p_child) = MPM%H(p_parent) + MPM%GradH(2, p_parent) * deta * (0.5 * dy)
  endif

  MPM%H(p_child) = max(MPM%H(p_child), MPM%H_min)

  if (MPM%length_update_type == L_CORNERS) then
    MPM%GVolume(p_child) = MPM%Lx(p_child) * MPM%Ly(p_child)
  else
   MPM%GVolume(p_child) = vol_frac * MPM%GVolume(p_parent)
  endif

  MPM%Fdef(:,p_child) = MPM%Fdef(:,p_parent)
  MPM%Bp(:,p_child) = MPM%Bp(:,p_parent)
  MPM%GradVel(:,p_child) = MPM%GradVel(:,p_parent)
  MPM%GradZs(:,p_child)  = MPM%GradZs(:,p_parent)
  MPM%GradH(:,p_child)  = MPM%GradH(:,p_parent)
  MPM%AGLen(p_child)            = MPM%AGLen(p_parent)
  MPM%eta_visc(p_child)         = MPM%eta_visc(p_parent)
  MPM%basal_trac_p(p_child)     = MPM%basal_trac_p(p_parent)
  MPM%newton_drag_coef_p(p_child) = MPM%newton_drag_coef_p(p_parent)
  MPM%up_mid(p_child)           = MPM%up_mid(p_parent)
  MPM%vp_mid(p_child)           = MPM%vp_mid(p_parent)
  MPM%smb(p_child)              = MPM%smb(p_parent)
  MPM%newton_visc_factor(p_child) = MPM%newton_visc_factor(p_parent)
  MPM%status(p_child) = ALIVE
  MPM%pid(p_child)    = MPM%pid_counter
  MPM%pid_counter = MPM%pid_counter + 1

end subroutine split_particle_child

! ============================================================
!> Re-assign particles to their owning cells after position updates.
!! Particles that have left the computational domain are marked DEAD.
!! Particles that have crossed PE boundaries are handled by simple
!! compaction here (full MPI migration is a future extension).
subroutine MPM_migrate(MPM, ISS, G)
  type(MPM_CS),          intent(inout) :: MPM  !< MPM control structure
  type(ice_shelf_state), intent(in)    :: ISS  !< Ice shelf state
  type(ocean_grid_type), intent(in)    :: G    !< Ocean grid

  integer :: p, ic_new, jc_new, ic_old, jc_old
  real :: xi_p, eta_p, xi_new, eta_new
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  ! Reset cell counts before re-assignment
  MPM%cell_count(:,:) = 0

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle

    ic_old = MPM%ci(p) ; jc_old = MPM%cj(p)
    xi_p  = MPM%xi(p)
    eta_p = MPM%eta(p)

    ! Determine new cell by tracking which cell the particle is in
    ic_new = ic_old ; jc_new = jc_old

    ! Walk in x
    do while (xi_p > 1.0 .and. ic_new <= ied)
      ic_new = ic_new + 1
      xi_p   = xi_p - 2.0
    enddo
    do while (xi_p < -1.0 .and. ic_new >= isd)
      ic_new = ic_new - 1
      xi_p   = xi_p + 2.0
    enddo

    ! Walk in y
    do while (eta_p > 1.0 .and. jc_new <= jed)
      jc_new = jc_new + 1
      eta_p  = eta_p - 2.0
    enddo
    do while (eta_p < -1.0 .and. jc_new >= jsd)
      jc_new = jc_new - 1
      eta_p  = eta_p + 2.0
    enddo

    ! Clamp updated local coords
    MPM%xi(p)  = xi_p
    MPM%eta(p) = eta_p

    ! Check if particle has left the domain or calved out
    if (ic_new < isc .or. ic_new > iec .or. &
        jc_new < jsc .or. jc_new > jec) then
      MPM%status(p) = DEAD
      cycle
    endif

    ! Check hmask: if cell is permanently ice-free (e.g. open ocean past calving mask)
    if (ISS%hmask(ic_new, jc_new) == 0 .and. ic_new /= ic_old) then
      ! Particle has entered an ice-free cell; mark for calving
      ! (front advance is permitted in hmask=0 cells adjacent to ice)
    endif

    MPM%ci(p) = ic_new ; MPM%cj(p) = jc_new
    MPM%cell_count(ic_new, jc_new) = MPM%cell_count(ic_new, jc_new) + 1
  enddo

  ! Compact DEAD particles (swap-with-last)
  call MPM_compact(MPM)

  ! Rebuild cell connectivity
  call MPM_build_cell_list(MPM, isd, ied, jsd, jed)

end subroutine MPM_migrate

! ============================================================
!> Remove DEAD particles by compaction (swap-with-last strategy).
subroutine MPM_compact(MPM)
  type(MPM_CS), intent(inout) :: MPM

  integer :: p, last

  last = MPM%n_active
  p = 1
  do while (p <= last)
    if (MPM%status(p) == DEAD) then
      ! Swap with last alive particle
      call swap_particles(MPM, p, last)
      last = last - 1
      ! Don't advance p - check the swapped-in particle next iteration
    else
      p = p + 1
    endif
  enddo
  MPM%n_active = last

end subroutine MPM_compact

! ============================================================
!> Enforce analytical Weertman thickness and velocity on particles in
!! hmask=3 cells for the 1D flow-band verification test (Huth et al. 2021).
!!
!! Extension domain (x <= inflow_x): H = H0, u = v0
!! First active cell (x > inflow_x, hmask=3): H = analytical, u = Q0/H
!!
!! The analytical steady-state profile is:
!!   alpha = A * (rho_i * g * (1 - rho_i/rho_w) / 4)^n
!!   H(x_rel) = (4*alpha*x_rel/Q0 + 1/H0^4)^(-0.25)
!!   u(x_rel) = Q0 / H(x_rel)
!! where x_rel = x - inflow_x and Q0 = H0 * v0.
subroutine MPM_enforce_1d_test(MPM, ISS, G)
  type(MPM_CS),          intent(inout) :: MPM   !< MPM control structure
  type(ice_shelf_state), intent(in)    :: ISS   !< Ice shelf state
  type(ocean_grid_type), intent(in)    :: G     !< Ocean grid

  integer :: p, ic, jc
  real :: x_km, x_rel_m, dx_km
  real :: H0, v0, Q0, A_glen, n_g
  real :: drive, alpha, m1, m2, H_si, u_si

  ! Hardcoded Huth et al. 2021 physical constants (SI)
  real, parameter :: rhoi_si = 910.0       ! [kg m-3]
  real, parameter :: rhow_si = 1028.0      ! [kg m-3]
  real, parameter :: grav_si = 9.81        ! [m s-2]

  if (.not. MPM%do_1d_test) return

  H0     = MPM%test1d_H0_si        ! [m]
  v0     = MPM%test1d_v0_si         ! [m s-1]
  Q0     = H0 * v0                  ! [m2 s-1]
  A_glen = MPM%AGlen_const          ! [Pa-n s-1] (raw SI)
  n_g    = MPM%n_glen               ! Glen exponent (typically 3)

  ! Weertman driving stress per unit thickness [Pa m-1]
  drive = rhoi_si * grav_si * (1.0 - rhoi_si / rhow_si) / 4.0
  ! alpha [m-n s-1]: strain-rate constant without H factor
  alpha = A_glen * drive**n_g
  ! Coefficients for H(x) = (m1*x + m2)^(-0.25) (valid for n=3)
  m1 = 4.0 * alpha / Q0   ! [m-4]
  m2 = 1.0 / (H0**4)      ! [m-4]

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    if (ISS%hmask(ic,jc) /= 1 .and. ISS%hmask(ic,jc) /= 3) cycle

    ! Physical x-position of particle [km]
    ! geoLonT is cell centre in km; xi is local coord in [-1,1]
    dx_km = G%dxT(ic,jc) * MPM%L_to_m * 1.0e-3   ! cell width in km
    x_km  = G%geoLonT(ic,jc) + MPM%xi(p) * 0.5 * dx_km

    ! Distance from inflow boundary [m]
    x_rel_m = (x_km - MPM%test1d_inflow_x_km) * 1.0e3

    if (x_rel_m <= 0.0) then
      ! Extension domain: enforce inflow values
      MPM%H(p)  = H0 * MPM%m_to_Z
      MPM%up(p) = v0 * MPM%m_s_to_L_T
      MPM%vp(p) = 0.0
    elseif (x_rel_m <= 5000.0) then
      ! First active cell: enforce analytical Weertman profile
      H_si = (m1 * x_rel_m + m2)**(-0.25)
      !u_si = Q0 / H_si
      MPM%H(p)  = H_si * MPM%m_to_Z
      !MPM%up(p) = u_si * MPM%m_s_to_L_T
      !MPM%vp(p) = 0.0
    endif
  enddo

end subroutine MPM_enforce_1d_test

! ============================================================
!> Write particle state to a per-PE VTK Unstructured Grid (.vtu) file for
!! visualisation in ParaView.  The root PE also writes a parallel collection
!! file (.pvtu) that references all per-PE files.
!!
!! Coordinates are output in metres (L ~> m at default unit scaling).
!! The file naming convention is:
!!   <output_dir>/mpm_<step>.pvtu        (root PE only)
!!   <output_dir>/mpm_<step>_pe<pe>.vtu  (every PE)
!!
!! @param MPM       MPM control structure
!! @param G         Ocean grid (provides dxT, dyT, idg_offset, jdg_offset)
!! @param step      Simulation step counter used in file names
!! @param output_dir Directory for output files (no trailing slash)
subroutine MPM_write_vtu(MPM, G, model_time_total, output_dir)
  type(MPM_CS),          intent(in) :: MPM              !< MPM control structure
  type(ocean_grid_type), intent(in) :: G                !< Ocean grid
  real,                  intent(in) :: model_time_total !< Cumulative model time [T ~> s]; used for filename
  character(len=*),      intent(in) :: output_dir       !< Directory for VTU files

  integer :: p, np, n_alive, iunit, step_num, k
  integer :: pe, npe
  real    :: xp, yp
  character(len=256) :: vtu_fname, pvtu_fname, pe_basename
  character(len=20)  :: step_str, pe_str
  integer :: ic, jc
  ! Convert model velocity [L T-1] to m year-1:
  !   vel [m a-1] = vel [L T-1] / m_s_to_L_T * seconds_per_year
  real, parameter :: seconds_per_year = 365.0 * 86400.0
  real :: to_m_a

  to_m_a = seconds_per_year / MPM%m_s_to_L_T
  pe  = PE_here()
  npe = num_PEs()

  ! Step number derived from model time — restart-safe, same file produced for same time
  if (MPM%vtu_dt > 0.0) then
    step_num = nint(model_time_total / MPM%vtu_dt)
  else
    step_num = 0
  endif

  ! Ensure the output directory exists (root PE creates it; benign if it already exists)
  if (is_root_pe()) call execute_command_line("mkdir -p " // trim(output_dir), wait=.true.)

  ! --- Build file names ---
  write(step_str, '(i0.6)') step_num
  write(pe_str,   '(i0.4)') pe

  vtu_fname  = trim(output_dir) // "/mpm_" // trim(step_str) // "_pe" // trim(pe_str) // ".vtu"
  pvtu_fname = trim(output_dir) // "/mpm_" // trim(step_str) // ".pvtu"

  ! Count ALIVE particles on this PE
  n_alive = 0
  do p = 1, MPM%n_active
    if (MPM%status(p) == ALIVE) n_alive = n_alive + 1
  enddo

  ! -------------------------------------------------------
  ! Write per-PE .vtu file (ASCII VTK XML, vertex elements)
  ! -------------------------------------------------------
  open(newunit=iunit, file=trim(vtu_fname), status='replace', action='write')

  write(iunit,'(a)') '<?xml version="1.0"?>'
  write(iunit,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(iunit,'(a)') '  <UnstructuredGrid>'
  write(iunit,'(a,i0,a,i0,a)') '    <Piece NumberOfPoints="', n_alive, &
                                '" NumberOfCells="', n_alive, '">'

  ! Points (x, y, z=0)
  write(iunit,'(a)') '      <Points>'
  write(iunit,'(a)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    xp = (real(ic + G%idg_offset) - 0.5) * G%dxT(ic,jc) + MPM%xi(p)  * G%dxT(ic,jc) * 0.5
    yp = (real(jc + G%jdg_offset) - 0.5) * G%dyT(ic,jc) + MPM%eta(p) * G%dyT(ic,jc) * 0.5
    write(iunit,'(3es20.10)') xp, yp, 0.0d0
  enddo
  write(iunit,'(a)') '        </DataArray>'
  write(iunit,'(a)') '      </Points>'

  ! Cells: one vertex cell per particle
  write(iunit,'(a)') '      <Cells>'
  write(iunit,'(a)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
  np = 0
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(i0)') np
    np = np + 1
  enddo
  write(iunit,'(a)') '        </DataArray>'
  write(iunit,'(a)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
  do p = 1, n_alive
    write(iunit,'(i0)') p
  enddo
  write(iunit,'(a)') '        </DataArray>'
  write(iunit,'(a)') '        <DataArray type="UInt8" Name="types" format="ascii">'
  do p = 1, n_alive
    write(iunit,'(i0)') 1   ! VTK_VERTEX
  enddo
  write(iunit,'(a)') '        </DataArray>'
  write(iunit,'(a)') '      </Cells>'

  ! Point data
  write(iunit,'(a)') '      <PointData>'

  ! H (thickness [Z ~> m])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="H_m" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%H(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! PVolume (current area [L2 ~> m2])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="PVolume_m2" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%PVolume(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! up (zonal particle velocity [m a-1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="up_ma" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%up(p) * to_m_a)
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! vp (meridional particle velocity [m a-1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="vp_ma" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%vp(p) * to_m_a)
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! up_g (pre-solve grid velocity, zonal [m a-1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="up_g_ma" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%up_g(p) * to_m_a)
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! vp_g (pre-solve grid velocity, meridional [m a-1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="vp_g_ma" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%vp_g(p) * to_m_a)
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! GVolume (reference area)
  write(iunit,'(a)') '        <DataArray type="Float32" Name="GVolume_m2" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%GVolume(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! Lx (current x-dimension)
  write(iunit,'(a)') '        <DataArray type="Float32" Name="Lx_m" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%Lx(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! Ly (current y-dimension)
  write(iunit,'(a)') '        <DataArray type="Float32" Name="Ly_m" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%Ly(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! Lx0 (reference x-dimension)
  write(iunit,'(a)') '        <DataArray type="Float32" Name="Lx0_m" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%Lx0(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! Ly0 (reference y-dimension)
  write(iunit,'(a)') '        <DataArray type="Float32" Name="Ly0_m" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%Ly0(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! strain_x
  write(iunit,'(a)') '        <DataArray type="Float32" Name="strain_x" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%strain_x(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! strain_y
  write(iunit,'(a)') '        <DataArray type="Float32" Name="strain_y" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%strain_y(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! Fdef components (4 scalar fields)
  do k = 1, 4
    write(iunit,'(a,i1,a)') '        <DataArray type="Float32" Name="Fdef', k, '" format="ascii">'
    do p = 1, MPM%n_active
      if (MPM%status(p) /= ALIVE) cycle
      write(iunit,'(es15.6)') real(MPM%Fdef(k,p))
    enddo
    write(iunit,'(a)') '        </DataArray>'
  enddo

  ! AGLen (rate factor)
  write(iunit,'(a)') '        <DataArray type="Float32" Name="AGLen" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%AGLen(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! smb (surface mass balance [Z T-1 ~> m s-1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="smb_ms" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%smb(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! eta_visc (effective viscosity [R L4 Z T-1 ~> kg m2 s-1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="eta_visc" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%eta_visc(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! ci (owning cell i-index)
  write(iunit,'(a)') '        <DataArray type="Int32" Name="ci" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(i0)') MPM%ci(p)
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! cj (owning cell j-index)
  write(iunit,'(a)') '        <DataArray type="Int32" Name="cj" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(i0)') MPM%cj(p)
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! xi (cell-local x coordinate [-1,1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="xi" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%xi(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! eta (cell-local y coordinate [-1,1])
  write(iunit,'(a)') '        <DataArray type="Float32" Name="eta" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(es15.6)') real(MPM%eta(p))
  enddo
  write(iunit,'(a)') '        </DataArray>'

  ! Bp (APIC affine matrix, 4 components)
  do k = 1, 4
    write(iunit,'(a,i1,a)') '        <DataArray type="Float32" Name="Bp', k, '" format="ascii">'
    do p = 1, MPM%n_active
      if (MPM%status(p) /= ALIVE) cycle
      write(iunit,'(es15.6)') real(MPM%Bp(k,p))
    enddo
    write(iunit,'(a)') '        </DataArray>'
  enddo

  ! pid (particle ID)
  write(iunit,'(a)') '        <DataArray type="Int64" Name="pid" format="ascii">'
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    write(iunit,'(i0)') MPM%pid(p)
  enddo
  write(iunit,'(a)') '        </DataArray>'

  write(iunit,'(a)') '      </PointData>'
  write(iunit,'(a)') '    </Piece>'
  write(iunit,'(a)') '  </UnstructuredGrid>'
  write(iunit,'(a)') '</VTKFile>'

  close(iunit)

  ! ---------------------------------------------------------
  ! Root PE writes a .pvtu parallel collection file
  ! ---------------------------------------------------------
  if (is_root_pe()) then
    open(newunit=iunit, file=trim(pvtu_fname), status='replace', action='write')

    write(iunit,'(a)') '<?xml version="1.0"?>'
    write(iunit,'(a)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iunit,'(a)') '  <PUnstructuredGrid GhostLevel="0">'
    write(iunit,'(a)') '    <PPoints>'
    write(iunit,'(a)') '      <PDataArray type="Float64" NumberOfComponents="3"/>'
    write(iunit,'(a)') '    </PPoints>'
    write(iunit,'(a)') '    <PPointData>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="H_m"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="PVolume_m2"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="up_ma"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="vp_ma"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="up_g_ma"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="vp_g_ma"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="GVolume_m2"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Lx_m"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Ly_m"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Lx0_m"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Ly0_m"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="strain_x"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="strain_y"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Fdef1"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Fdef2"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Fdef3"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Fdef4"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="AGLen"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="smb_ms"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="eta_visc"/>'
    write(iunit,'(a)') '      <PDataArray type="Int32"   Name="ci"/>'
    write(iunit,'(a)') '      <PDataArray type="Int32"   Name="cj"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="xi"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="eta"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Bp1"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Bp2"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Bp3"/>'
    write(iunit,'(a)') '      <PDataArray type="Float32" Name="Bp4"/>'
    write(iunit,'(a)') '      <PDataArray type="Int64"   Name="pid"/>'
    write(iunit,'(a)') '    </PPointData>'

    do pe = 0, npe - 1
      write(pe_str, '(i0.4)') pe
      pe_basename = "mpm_" // trim(step_str) // "_pe" // trim(pe_str) // ".vtu"
      write(iunit,'(a)') '    <Piece Source="' // trim(pe_basename) // '"/>'
    enddo

    write(iunit,'(a)') '  </PUnstructuredGrid>'
    write(iunit,'(a)') '</VTKFile>'

    close(iunit)
  endif

end subroutine MPM_write_vtu

! ============================================================
!> Save all ALIVE particle state to a per-PE binary restart file.
!!
!! Files are written to <restart_dir>/mpm_peXXXX.res (one per PE).
!! The format is a Fortran unformatted sequential file:
!!   Record 1 : integer(4) n_alive, integer(8) pid_counter
!!   Record 2 : real(8) xi(1:n_alive)
!!   Record 3 : real(8) eta(1:n_alive)
!!   Record 4 : integer(4) ci(1:n_alive)
!!   Record 5 : integer(4) cj(1:n_alive)
!!   Record 6 : real(8) up(1:n_alive)
!!   Record 7 : real(8) vp(1:n_alive)
!!   Record 8 : real(8) up_g(1:n_alive)
!!   Record 9 : real(8) vp_g(1:n_alive)
!!   Record 10: real(8) H(1:n_alive)
!!   Record 11: real(8) PVolume(1:n_alive)
!!   Record 12: real(8) GVolume(1:n_alive)
!!   Record 13: real(8) Lx(1:n_alive)
!!   Record 14: real(8) Ly(1:n_alive)
!!   Record 15: real(8) Lx0(1:n_alive)
!!   Record 16: real(8) Ly0(1:n_alive)
!!   Record 17: real(8) strain_x(1:n_alive)
!!   Record 18: real(8) strain_y(1:n_alive)
!!   Record 19: real(8) Fdef(4,1:n_alive)
!!   Record 20: real(8) AGLen(1:n_alive)
!!   Record 21: real(8) smb(1:n_alive)
!!   Record 22: integer(8) pid(1:n_alive)
subroutine MPM_save_restart(MPM, restart_dir)
  type(MPM_CS),     intent(in) :: MPM         !< MPM control structure
  character(len=*), intent(in) :: restart_dir !< Directory for restart files

  integer :: p, n_alive, iunit, ios
  character(len=256) :: fname
  character(len=20)  :: pe_str

  ! --- Count ALIVE particles ---
  n_alive = 0
  do p = 1, MPM%n_active
    if (MPM%status(p) == ALIVE) n_alive = n_alive + 1
  enddo

  ! --- Build filename ---
  write(pe_str, '(i0.4)') PE_here()
  fname = trim(restart_dir) // "/mpm_pe" // trim(pe_str) // ".res"

  ! --- Pack ALIVE particle arrays into temporary 1-D arrays ---
  ! To avoid holding large temporaries we write field-by-field using
  ! a compaction index array built inline.
  block
    integer,          allocatable :: idx(:)
    real,             allocatable :: rbuf(:)
    integer,          allocatable :: ibuf(:)
    integer(kind=8),  allocatable :: i8buf(:)
    real,             allocatable :: rbuf2(:,:)
    integer :: k, cnt

    allocate(idx(n_alive))
    cnt = 0
    do p = 1, MPM%n_active
      if (MPM%status(p) /= ALIVE) cycle
      cnt = cnt + 1
      idx(cnt) = p
    enddo

    allocate(rbuf(n_alive), ibuf(n_alive), i8buf(n_alive), rbuf2(4, n_alive))

    open(newunit=iunit, file=trim(fname), status='replace', action='write', &
         form='unformatted', iostat=ios)
    if (ios /= 0) call MOM_error(FATAL, "MPM_save_restart: cannot open " // trim(fname))

    ! Record 1: counts and VTU state
    write(iunit) n_alive, MPM%pid_counter, MPM%model_time_total, MPM%vtu_time_accum

    ! xi
    do k = 1, n_alive ; rbuf(k) = MPM%xi(idx(k))      ; enddo ; write(iunit) rbuf
    ! eta
    do k = 1, n_alive ; rbuf(k) = MPM%eta(idx(k))     ; enddo ; write(iunit) rbuf
    ! ci
    do k = 1, n_alive ; ibuf(k) = MPM%ci(idx(k))      ; enddo ; write(iunit) ibuf
    ! cj
    do k = 1, n_alive ; ibuf(k) = MPM%cj(idx(k))      ; enddo ; write(iunit) ibuf
    ! up
    do k = 1, n_alive ; rbuf(k) = MPM%up(idx(k))      ; enddo ; write(iunit) rbuf
    ! vp
    do k = 1, n_alive ; rbuf(k) = MPM%vp(idx(k))      ; enddo ; write(iunit) rbuf
    ! up_g
    do k = 1, n_alive ; rbuf(k) = MPM%up_g(idx(k))    ; enddo ; write(iunit) rbuf
    ! vp_g
    do k = 1, n_alive ; rbuf(k) = MPM%vp_g(idx(k))    ; enddo ; write(iunit) rbuf
    ! H
    do k = 1, n_alive ; rbuf(k) = MPM%H(idx(k))       ; enddo ; write(iunit) rbuf
    ! PVolume
    do k = 1, n_alive ; rbuf(k) = MPM%PVolume(idx(k)) ; enddo ; write(iunit) rbuf
    ! GVolume
    do k = 1, n_alive ; rbuf(k) = MPM%GVolume(idx(k)) ; enddo ; write(iunit) rbuf
    ! Lx
    do k = 1, n_alive ; rbuf(k) = MPM%Lx(idx(k))      ; enddo ; write(iunit) rbuf
    ! Ly
    do k = 1, n_alive ; rbuf(k) = MPM%Ly(idx(k))      ; enddo ; write(iunit) rbuf
    ! Lx0
    do k = 1, n_alive ; rbuf(k) = MPM%Lx0(idx(k))     ; enddo ; write(iunit) rbuf
    ! Ly0
    do k = 1, n_alive ; rbuf(k) = MPM%Ly0(idx(k))     ; enddo ; write(iunit) rbuf
    ! strain_x
    do k = 1, n_alive ; rbuf(k) = MPM%strain_x(idx(k)); enddo ; write(iunit) rbuf
    ! strain_y
    do k = 1, n_alive ; rbuf(k) = MPM%strain_y(idx(k)); enddo ; write(iunit) rbuf
    ! Fdef (4 components)
    do k = 1, n_alive ; rbuf2(:,k) = MPM%Fdef(:,idx(k)) ; enddo ; write(iunit) rbuf2
    ! AGLen
    do k = 1, n_alive ; rbuf(k) = MPM%AGLen(idx(k))   ; enddo ; write(iunit) rbuf
    ! smb
    do k = 1, n_alive ; rbuf(k) = MPM%smb(idx(k))     ; enddo ; write(iunit) rbuf
    ! pid
    do k = 1, n_alive ; i8buf(k) = MPM%pid(idx(k))    ; enddo ; write(iunit) i8buf

    close(iunit)
    deallocate(idx, rbuf, ibuf, i8buf, rbuf2)
  end block

  call MOM_mesg("MPM_save_restart: wrote " // trim(fname) // &
                " (" // trim(int_to_str(n_alive)) // " particles)", 3)

end subroutine MPM_save_restart

! ============================================================
!> Restore particle state from a per-PE binary restart file written by
!! MPM_save_restart.  Existing particle arrays are discarded; the MPM
!! state is rebuilt from the file and the cell-particle lists are
!! reconstructed.
subroutine MPM_restore_restart(MPM, G, restart_dir)
  type(MPM_CS),          intent(inout) :: MPM         !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G           !< Ocean grid (for bounds only)
  character(len=*),      intent(in)    :: restart_dir !< Directory for restart files

  integer :: n_alive, iunit, ios
  integer(kind=8) :: pid_ctr
  character(len=256) :: fname
  character(len=20)  :: pe_str

  write(pe_str, '(i0.4)') PE_here()
  fname = trim(restart_dir) // "/mpm_pe" // trim(pe_str) // ".res"

  open(newunit=iunit, file=trim(fname), status='old', action='read', &
       form='unformatted', iostat=ios)
  if (ios /= 0) call MOM_error(FATAL, &
       "MPM_restore_restart: cannot open " // trim(fname))

  ! Read record 1: n_alive, pid_counter, and VTU state
  read(iunit) n_alive, pid_ctr, MPM%model_time_total, MPM%vtu_time_accum

  ! Grow arrays until n_alloc >= n_alive
  do while (n_alive > MPM%n_alloc)
    call MPM_grow_arrays(MPM, G%isd, G%ied, G%jsd, G%jed, &
                              G%IsdB, G%IedB, G%JsdB, G%JedB)
  enddo
  MPM%n_active    = n_alive
  MPM%pid_counter = pid_ctr

  ! Read particle arrays directly (all particles are ALIVE after restart)
  read(iunit) MPM%xi(1:n_alive)
  read(iunit) MPM%eta(1:n_alive)
  read(iunit) MPM%ci(1:n_alive)
  read(iunit) MPM%cj(1:n_alive)
  read(iunit) MPM%up(1:n_alive)
  read(iunit) MPM%vp(1:n_alive)
  read(iunit) MPM%up_g(1:n_alive)
  read(iunit) MPM%vp_g(1:n_alive)
  read(iunit) MPM%H(1:n_alive)
  read(iunit) MPM%PVolume(1:n_alive)
  read(iunit) MPM%GVolume(1:n_alive)
  read(iunit) MPM%Lx(1:n_alive)
  read(iunit) MPM%Ly(1:n_alive)
  read(iunit) MPM%Lx0(1:n_alive)
  read(iunit) MPM%Ly0(1:n_alive)
  read(iunit) MPM%strain_x(1:n_alive)
  read(iunit) MPM%strain_y(1:n_alive)
  read(iunit) MPM%Fdef(:,1:n_alive)
  read(iunit) MPM%AGLen(1:n_alive)
  read(iunit) MPM%smb(1:n_alive)
  read(iunit) MPM%pid(1:n_alive)

  close(iunit)

  ! Mark all restored particles as ALIVE and rebuild cell counts from ci/cj
  MPM%status(1:n_alive) = ALIVE
  MPM%cell_count(:,:) = 0
  block
    integer :: p, ic, jc
    do p = 1, n_alive
      ic = MPM%ci(p) ; jc = MPM%cj(p)
      MPM%cell_count(ic, jc) = MPM%cell_count(ic, jc) + 1
    enddo
  end block

  ! Rebuild cell-particle connectivity
  call MPM_build_cell_list(MPM, G%isd, G%ied, G%jsd, G%jed)

  call MOM_mesg("MPM_restore_restart: read " // trim(fname) // &
                " (" // trim(int_to_str(n_alive)) // " particles)", 3)

end subroutine MPM_restore_restart

! ============================================================
!> Deallocate all MPM arrays.
subroutine MPM_end(MPM)
  type(MPM_CS), intent(inout) :: MPM

  if (allocated(MPM%xi)) deallocate(MPM%xi)
  if (allocated(MPM%eta)) deallocate(MPM%eta)
  if (allocated(MPM%ci)) deallocate(MPM%ci)
  if (allocated(MPM%cj)) deallocate(MPM%cj)
  if (allocated(MPM%up)) deallocate(MPM%up)
  if (allocated(MPM%vp)) deallocate(MPM%vp)
  if (allocated(MPM%up_g)) deallocate(MPM%up_g)
  if (allocated(MPM%vp_g)) deallocate(MPM%vp_g)
  if (allocated(MPM%H)) deallocate(MPM%H)
  if (allocated(MPM%PVolume)) deallocate(MPM%PVolume)
  if (allocated(MPM%GVolume)) deallocate(MPM%GVolume)
  if (allocated(MPM%Lx)) deallocate(MPM%Lx)
  if (allocated(MPM%Ly)) deallocate(MPM%Ly)
  if (allocated(MPM%Lx0)) deallocate(MPM%Lx0)
  if (allocated(MPM%Ly0)) deallocate(MPM%Ly0)
  if (allocated(MPM%strain_x)) deallocate(MPM%strain_x)
  if (allocated(MPM%strain_y)) deallocate(MPM%strain_y)
  if (allocated(MPM%Fdef)) deallocate(MPM%Fdef)
  if (allocated(MPM%GradVel)) deallocate(MPM%GradVel)
  if (allocated(MPM%GradZs)) deallocate(MPM%GradZs)
  if (allocated(MPM%GradH))  deallocate(MPM%GradH)
  if (allocated(MPM%AGLen)) deallocate(MPM%AGLen)
  if (allocated(MPM%eta_visc)) deallocate(MPM%eta_visc)
  if (allocated(MPM%newton_visc_factor)) deallocate(MPM%newton_visc_factor)
  if (allocated(MPM%basal_trac_p)) deallocate(MPM%basal_trac_p)
  if (allocated(MPM%newton_drag_coef_p)) deallocate(MPM%newton_drag_coef_p)
  if (allocated(MPM%up_mid)) deallocate(MPM%up_mid)
  if (allocated(MPM%vp_mid)) deallocate(MPM%vp_mid)
  if (allocated(MPM%Bp)) deallocate(MPM%Bp)
  if (allocated(MPM%smb)) deallocate(MPM%smb)
  if (allocated(MPM%is_front_cell)) deallocate(MPM%is_front_cell)
  if (allocated(MPM%status)) deallocate(MPM%status)
  if (allocated(MPM%pid)) deallocate(MPM%pid)
  if (allocated(MPM%cell_list)) deallocate(MPM%cell_list)
  if (allocated(MPM%cell_count)) deallocate(MPM%cell_count)
  if (allocated(MPM%cell_start)) deallocate(MPM%cell_start)
  if (allocated(MPM%H_node)) deallocate(MPM%H_node)
  if (allocated(MPM%Zs_node)) deallocate(MPM%Zs_node)
  if (allocated(MPM%bed_node)) deallocate(MPM%bed_node)
  if (allocated(MPM%h_node_mask)) deallocate(MPM%h_node_mask)
  if (allocated(MPM%h_node_bdry_val)) deallocate(MPM%h_node_bdry_val)
  if (allocated(MPM%area_shelf_h)) deallocate(MPM%area_shelf_h)
  if (allocated(MPM%reweight)) deallocate(MPM%reweight)

  MPM%initialized = .false.

end subroutine MPM_end

! ============================================================
!> Compute per-particle effective viscosity from current velocity gradients.
!! Called each Picard iteration. Uses Glen's flow law:
!!   eta(p) = 0.5 * H(p) * A^(-1/n) * eps_e^((1-n)/n)
!! where eps_e is the effective strain rate invariant.
subroutine MPM_compute_visc(MPM, G, US, u_shelf, v_shelf)
  type(MPM_CS),          intent(inout) :: MPM
  type(ocean_grid_type), intent(in)    :: G
  type(unit_scale_type), intent(in)    :: US  !< Unit scaling factors
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), intent(in) :: u_shelf, v_shelf

  integer :: p, ic, jc, k_n, n_gimp, I_Bu, J_Bu, IsdB, IedB, JsdB, JedB
  real :: xi_p, eta_p, dx, dy
  real :: Ngimp(9), dNdx_gimp(9), dNdy_gimp(9)
  integer :: di_gimp(9), dj_gimp(9)
  real :: ux, vy, uy, vx_g
  real :: eps_e2, eps_e, visc_coef, n_g, one_over_n, eps_min_sq

  n_g = MPM%n_glen
  one_over_n = 1.0 / n_g
  eps_min_sq = MPM%eps_glen_min * MPM%eps_glen_min
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    dx = G%dxT(ic,jc) ; dy = G%dyT(ic,jc)
    xi_p = MPM%xi(p) ; eta_p = MPM%eta(p)

    if (MPM%basis_type == BASIS_GIMP) then
      call gimp_nodes(xi_p, eta_p, 0.5*MPM%Lx(p), 0.5*MPM%Ly(p), dx, dy, &
                      Ngimp, dNdx_gimp, dNdy_gimp, di_gimp, dj_gimp, n_gimp)
    else
      call smpm_nodes(xi_p, eta_p, dx, dy, &
                      Ngimp, dNdx_gimp, dNdy_gimp, di_gimp, dj_gimp, n_gimp)
    endif

    ux = 0.0 ; vy = 0.0 ; uy = 0.0 ; vx_g = 0.0
    do k_n = 1, n_gimp
      I_Bu = ic - 1 + di_gimp(k_n)
      J_Bu = jc - 1 + dj_gimp(k_n)
      if (I_Bu < IsdB .or. I_Bu > IedB) cycle
      if (J_Bu < JsdB .or. J_Bu > JedB) cycle
      ux   = ux   + dNdx_gimp(k_n) * u_shelf(I_Bu, J_Bu)
      vy   = vy   + dNdy_gimp(k_n) * v_shelf(I_Bu, J_Bu)
      uy   = uy   + dNdy_gimp(k_n) * u_shelf(I_Bu, J_Bu)
      vx_g = vx_g + dNdx_gimp(k_n) * v_shelf(I_Bu, J_Bu)
    enddo

    ! Update stored velocity gradients
    MPM%GradVel(1,p) = ux
    MPM%GradVel(2,p) = vy
    MPM%GradVel(3,p) = uy
    MPM%GradVel(4,p) = vx_g

    ! Effective strain rate: eps_e^2 = ux^2 + vy^2 + ux*vy + 0.25*(uy+vx)^2 + eps_min^2
    ! The eps_min^2 regularization matches calc_shelf_visc (line 4799)
    eps_e2 = ux*ux + vy*vy + ux*vy + 0.25*(uy + vx_g)*(uy + vx_g) &
           + MPM%eps_glen_min*MPM%eps_glen_min

    ! Visc_coef = A^(-1/n); AGLen(p) stores the raw A [Pa^-n s^-1]
    visc_coef = MPM%AGLen(p)**(-one_over_n)

    ! eta = 0.5 * H * visc_coef * eps_e2^((1-n)/(2n)) with unit conversion.
    ! This matches calc_shelf_visc but WITHOUT area scaling (CG_action_MPM
    ! applies the area weight via PVolume*reweight instead).
    ! eps_e2 is in [T^-2]; convert to [s^-2] with US%s_to_T^2 before the
    ! fractional power, then convert [Pa s] result to [R L^2 T^-1].
    MPM%eta_visc(p) = 0.5 * MPM%H(p) * visc_coef * &
        (US%s_to_T**2 * eps_e2)**((1.0 - n_g) / (2.0 * n_g)) * &
        (US%Pa_to_RL2_T2 * US%s_to_T)

    ! Newton viscosity correction factor (same formula as calc_shelf_visc):
    !   nvf = (0.5*(1/n - 1)) / eps_e2 * eta_visc
    ! Uses eps_e2 WITHOUT the eps_min^2 floor (same convention as FEM):
    eps_e2 = ux*ux + vy*vy + ux*vy + 0.25*(uy + vx_g)*(uy + vx_g)  ! no floor
    MPM%newton_visc_factor(p) = (0.5 * (one_over_n - 1.0) / &
        (eps_e2 + eps_min_sq)) * MPM%eta_visc(p)
  enddo

end subroutine MPM_compute_visc

! ============================================================
!> Compute is_front_cell flag: a cell is a "front cell" if hmask==1 and
!! any of its 4 face-neighbors has hmask==0. Inflow BC cells (hmask==3)
!! are never classified as front cells.
subroutine MPM_compute_is_front_cell(MPM, ISS, G, u_node_mask)
  type(MPM_CS),          intent(inout) :: MPM
  type(ice_shelf_state), intent(in)    :: ISS
  type(ocean_grid_type), intent(in)    :: G
  real, dimension(SZDIB_(G),SZDJB_(G)), intent(in) :: u_node_mask !< MPM u-velocity Dirichlet mask
                                              !! at Bu corners; 1=prescribed, 0=free [nondim]

  integer :: i, j, isc, iec, jsc, jec

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec

  MPM%is_front_cell(:,:) = .false.

  do j = jsc, jec ; do i = isc, iec
    if (ISS%hmask(i,j) /= 1.0) cycle
    ! Cells whose 4 corners are all Dirichlet velocity nodes are inflow BC cells:
    ! they do not adjoin a calving front and must not be classified as front cells.
    if (u_node_mask(i-1,j-1) == 1.0 .and. u_node_mask(i,j-1) == 1.0 .and. &
        u_node_mask(i-1,j  ) == 1.0 .and. u_node_mask(i,j  ) == 1.0) cycle
    ! Check all 4 face neighbours.  Halo cells at the global domain boundary
    ! have hmask=0 (ocean) after pass_var, so a cell at the calving-front edge
    ! is correctly identified as a front cell by its eastern halo neighbour.
    if (ISS%hmask(i-1,j) == 0.0 .or. &
        ISS%hmask(i+1,j) == 0.0 .or. &
        ISS%hmask(i,j-1) == 0.0 .or. &
        ISS%hmask(i,j+1) == 0.0) then
      MPM%is_front_cell(i,j) = .true.
    endif
  enddo ; enddo

end subroutine MPM_compute_is_front_cell

! ============================================================
!> Store pre-solve grid velocity at each particle for FLIP update.
!! Must be called before the SSA solve so that G2P after the solve
!! can compute the velocity increment.
subroutine MPM_store_pre_solve_vel(MPM, G, u_shelf, v_shelf)
  type(MPM_CS),          intent(inout) :: MPM
  type(ocean_grid_type), intent(in)    :: G
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), intent(in) :: u_shelf, v_shelf

  integer :: p, ic, jc, k_n, n_gimp, I_Bu, J_Bu, IsdB, IedB, JsdB, JedB
  real :: N_g(9), dNdx_g(9), dNdy_g(9)
  integer :: di_g(9), dj_g(9)
  real :: u_g, v_g, dx, dy

  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    dx = G%dxT(ic,jc) ; dy = G%dyT(ic,jc)
    if (MPM%basis_type == BASIS_GIMP) then
      call gimp_nodes(MPM%xi(p), MPM%eta(p), 0.5*MPM%Lx(p), 0.5*MPM%Ly(p), dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    else
      call smpm_nodes(MPM%xi(p), MPM%eta(p), dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    endif

    u_g = 0.0 ; v_g = 0.0
    do k_n = 1, n_gimp
      I_Bu = ic - 1 + di_g(k_n)
      J_Bu = jc - 1 + dj_g(k_n)
      if (I_Bu < IsdB .or. I_Bu > IedB) cycle
      if (J_Bu < JsdB .or. J_Bu > JedB) cycle
      u_g = u_g + N_g(k_n) * u_shelf(I_Bu, J_Bu)
      v_g = v_g + N_g(k_n) * v_shelf(I_Bu, J_Bu)
    enddo

    MPM%up_g(p) = u_g
    MPM%vp_g(p) = v_g
  enddo

end subroutine MPM_store_pre_solve_vel

! ============================================================
!> Momentum-conserving P2G velocity mapping before SSA solve (Huth 2021, Eq 27).
!!
!! Computes: v_I = (Σ_p m_p * N_Ip * v_p_eff) / (Σ_p m_p * N_Ip)
!! where m_p = rho_ice * H_p * PVolume_p and, for APIC:
!!   v_p_eff = v_p + C_p * (x_I - x_p)
!!   C_p = B_p * diag(4/dx², 4/dy²)  (for bilinear shape functions)
!!
!! Only overwrites free nodes (umask==1). Dirichlet nodes (ufacemask==5)
!! already have correct velocity from boundary conditions.
!!
!! Must be called BEFORE MPM_store_pre_solve_vel so that the momentum-mapped
!! velocity becomes v_I^m in the FLIP update.
subroutine MPM_P2G_velocity(MPM, G, u_shlf, v_shlf, umask, vmask, rho_ice)
  type(MPM_CS),          intent(inout) :: MPM     !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G       !< Ocean grid
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(inout) :: u_shlf  !< Ice velocity updated at free nodes
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(inout) :: v_shlf  !< Ice velocity updated at free nodes
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(in)    :: umask   !< u-velocity mask (1=free, 5=Dirichlet, -2=ocean)
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(in)    :: vmask   !< v-velocity mask
  real,                  intent(in)    :: rho_ice !< Ice density [R ~> kg m-3]

  integer :: p, ic, jc, k_n, n_gimp, I_Bu, J_Bu, IsdB, IedB, JsdB, JedB, I, J
  real :: N_g(9), dNdx_g(9), dNdy_g(9)
  integer :: di_g(9), dj_g(9)
  real :: m_p, xi_p, eta_p, dx, dy, dx_half, dy_half
  real :: u_eff, v_eff  ! APIC-augmented particle velocity for P2G
  real :: dx_I_k, dy_I_k  ! Node-to-particle offset [L]
  real :: inv_dx2, inv_dy2  ! 4/dx², 4/dy²
  real, allocatable :: num_u(:,:), den_u(:,:)  ! Numerator/denominator accumulators
  real, allocatable :: num_v(:,:), den_v(:,:)

  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  allocate(num_u(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(den_u(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(num_v(IsdB:IedB, JsdB:JedB), source=0.0)
  allocate(den_v(IsdB:IedB, JsdB:JedB), source=0.0)

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    xi_p  = MPM%xi(p)
    eta_p = MPM%eta(p)
    dx = G%dxT(ic,jc) ; dy = G%dyT(ic,jc)
    dx_half = 0.5 * dx ; dy_half = 0.5 * dy

    if (MPM%basis_type == BASIS_GIMP) then
      call gimp_nodes(xi_p, eta_p, 0.5*MPM%Lx(p), 0.5*MPM%Ly(p), dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    else
      call smpm_nodes(xi_p, eta_p, dx, dy, &
                      N_g, dNdx_g, dNdy_g, di_g, dj_g, n_gimp)
    endif

    m_p = rho_ice * MPM%H(p) * MPM%PVolume(p)

    if (MPM%transfer_type == TRANSFER_APIC) then
      ! APIC P2G: v_eff = v_p + C_p * (x_I - x_p),  C_p = B_p * diag(4/dx², 4/dy²)
      inv_dx2 = 4.0 / (dx * dx)
      inv_dy2 = 4.0 / (dy * dy)
      do k_n = 1, n_gimp
        I_Bu = ic - 1 + di_g(k_n)
        J_Bu = jc - 1 + dj_g(k_n)
        if (I_Bu < IsdB .or. I_Bu > IedB) cycle
        if (J_Bu < JsdB .or. J_Bu > JedB) cycle
        ! Node-to-particle offset: x_I - x_p = (xi_I - xi_p)*dx_half
        ! where xi_I = 2*di_g - 1 (= -1,+1,+3 for di_g = 0,1,2)
        dx_I_k = (-1.0 + 2.0*real(di_g(k_n)) - xi_p)  * dx_half
        dy_I_k = (-1.0 + 2.0*real(dj_g(k_n)) - eta_p) * dy_half
        u_eff = MPM%up(p) + (MPM%Bp(1,p)*dx_I_k*inv_dx2 + MPM%Bp(2,p)*dy_I_k*inv_dy2)
        v_eff = MPM%vp(p) + (MPM%Bp(3,p)*dx_I_k*inv_dx2 + MPM%Bp(4,p)*dy_I_k*inv_dy2)
        num_u(I_Bu, J_Bu) = num_u(I_Bu, J_Bu) + m_p * N_g(k_n) * u_eff
        num_v(I_Bu, J_Bu) = num_v(I_Bu, J_Bu) + m_p * N_g(k_n) * v_eff
        den_u(I_Bu, J_Bu) = den_u(I_Bu, J_Bu) + m_p * N_g(k_n)
        den_v(I_Bu, J_Bu) = den_v(I_Bu, J_Bu) + m_p * N_g(k_n)
      enddo
    else
      ! FLIP/PIC: standard Eq 27 — constant velocity per particle
      do k_n = 1, n_gimp
        I_Bu = ic - 1 + di_g(k_n)
        J_Bu = jc - 1 + dj_g(k_n)
        if (I_Bu < IsdB .or. I_Bu > IedB) cycle
        if (J_Bu < JsdB .or. J_Bu > JedB) cycle
        num_u(I_Bu, J_Bu) = num_u(I_Bu, J_Bu) + m_p * N_g(k_n) * MPM%up(p)
        num_v(I_Bu, J_Bu) = num_v(I_Bu, J_Bu) + m_p * N_g(k_n) * MPM%vp(p)
        den_u(I_Bu, J_Bu) = den_u(I_Bu, J_Bu) + m_p * N_g(k_n)
        den_v(I_Bu, J_Bu) = den_v(I_Bu, J_Bu) + m_p * N_g(k_n)
      enddo
    endif
  enddo

  ! Halo exchange so corner nodes on PE boundaries receive contributions from neighbours
  call pass_var(num_u, G%domain, position=CORNER)
  call pass_var(den_u, G%domain, position=CORNER)
  call pass_var(num_v, G%domain, position=CORNER)
  call pass_var(den_v, G%domain, position=CORNER)

  ! Normalize and overwrite free nodes only
  ! umask/vmask at Bu corners: 1.0 = free SSA node, 5.0 = Dirichlet, -2.0 = ocean
  ! Only update where umask==1 to preserve Dirichlet BCs?
  do J = JsdB, JedB ; do I = IsdB, IedB
    if (umask(I,J) >= 1.0 .and. den_u(I,J) > 0.0) &
      u_shlf(I,J) = num_u(I,J) / den_u(I,J)
    if (vmask(I,J) >= 1.0 .and. den_v(I,J) > 0.0) &
      v_shlf(I,J) = num_v(I,J) / den_v(I,J)
  enddo ; enddo

  deallocate(num_u, den_u, num_v, den_v)

end subroutine MPM_P2G_velocity

! ============================================================
!> Compute per-particle linearized basal drag coefficient.
!! Called once per outer Picard iteration from ice_shelf_solve_outer,
!! analogous to calc_shelf_taub for FEM cells.
!!
!! For floating particles: basal_trac_p = 0.
!! For grounded particles: computes a linearized beta coefficient
!! using the same power-law friction law as calc_shelf_taub.
!! Both 1D test cases (MPM_weertman, MPM_front_advance_v2) are fully floating,
!! so this routine will set basal_trac_p = 0 for all particles.
subroutine MPM_compute_basal_trac(MPM, G, US, bathyT, u_shlf, v_shlf)
  type(MPM_CS),          intent(inout) :: MPM     !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G       !< Ocean grid
  type(unit_scale_type), intent(in)    :: US      !< Unit scaling
  real, dimension(G%isd:G%ied, G%jsd:G%jed), intent(in) :: bathyT  !< Bed depth [Z ~> m] (positive down)
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), intent(in) :: u_shlf !< Ice velocity [L T-1]
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), intent(in) :: v_shlf !< Ice velocity [L T-1]

  integer :: p, ic, jc
  real :: draft         !< Ice draft [Z ~> m]
  real :: N(4)          !< Shape function values [nondim]
  real :: up_mid_p      !< Midpoint u-velocity at particle [L T-1 ~> m s-1]
  real :: vp_mid_p      !< Midpoint v-velocity at particle [L T-1 ~> m s-1]

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)

    ! Interpolate current Picard iterate velocity to particle position.
    ! This is stored as the midpoint velocity for the Newton basal drag correction,
    ! matching what calc_shelf_taub does with newton_umid/vmid for the FEM cells.
    call smpm_shape(MPM%xi(p), MPM%eta(p), N)
    up_mid_p = N(1)*u_shlf(ic-1,jc-1) + N(2)*u_shlf(ic,jc-1) + &
               N(3)*u_shlf(ic-1,jc  ) + N(4)*u_shlf(ic,jc  )
    vp_mid_p = N(1)*v_shlf(ic-1,jc-1) + N(2)*v_shlf(ic,jc-1) + &
               N(3)*v_shlf(ic-1,jc  ) + N(4)*v_shlf(ic,jc  )
    MPM%up_mid(p) = up_mid_p
    MPM%vp_mid(p) = vp_mid_p

    ! Grounding check: ice is grounded if flotation draft > bed depth
    draft = (MPM%density_ice / MPM%density_ocean) * MPM%H(p)
    if (draft > bathyT(ic,jc)) then
      ! Grounded particle: stub — power-law basal drag not yet implemented.
      ! Both basal_trac_p and newton_drag_coef_p remain zero until a full
      ! Weertman/Coulomb sliding law is added here (matching calc_shelf_taub).
      MPM%basal_trac_p(p)       = 0.0
      MPM%newton_drag_coef_p(p) = 0.0
    else
      ! Floating: no basal drag
      MPM%basal_trac_p(p)       = 0.0
      MPM%newton_drag_coef_p(p) = 0.0
    endif
  enddo

end subroutine MPM_compute_basal_trac

! ============================================================
! Private helper: sMPM bilinear shape functions at (xi, eta) in [-1,1]²
! Corner ordering: 1=SW, 2=SE, 3=NW, 4=NE  (same as MOM6 Phi convention)
subroutine smpm_shape(xi, eta, N)
  real, intent(in)  :: xi, eta   !< Local coordinates [nondim]
  real, intent(out) :: N(4)      !< Shape function values [nondim]

  N(1) = 0.25 * (1.0 - xi) * (1.0 - eta)  ! SW
  N(2) = 0.25 * (1.0 + xi) * (1.0 - eta)  ! SE
  N(3) = 0.25 * (1.0 - xi) * (1.0 + eta)  ! NW
  N(4) = 0.25 * (1.0 + xi) * (1.0 + eta)  ! NE

end subroutine smpm_shape

! ============================================================
! Private helper: sMPM shape function gradients at (xi, eta).
! dNdx(k) = ∂N_k/∂x,  dNdy(k) = ∂N_k/∂y
! For a uniform Cartesian cell of physical size dx × dy.
subroutine smpm_grad(xi, eta, dx, dy, dNdx, dNdy)
  real, intent(in)  :: xi, eta     !< Local coordinates [nondim]
  real, intent(in)  :: dx, dy      !< Physical cell dimensions [L]
  real, intent(out) :: dNdx(4)     !< ∂N_k/∂x [L-1]
  real, intent(out) :: dNdy(4)     !< ∂N_k/∂y [L-1]

  real :: inv2dx, inv2dy

  inv2dx = 1.0 / (2.0 * dx)   ! = 1/(2*dx)
  inv2dy = 1.0 / (2.0 * dy)

  ! ∂N_k/∂x = (∂N_k/∂ξ) * (2/dx)
  dNdx(1) = -(1.0 - eta) * inv2dx  ! SW
  dNdx(2) = +(1.0 - eta) * inv2dx  ! SE
  dNdx(3) = -(1.0 + eta) * inv2dx  ! NW
  dNdx(4) = +(1.0 + eta) * inv2dx  ! NE

  ! ∂N_k/∂y = (∂N_k/∂η) * (2/dy)
  dNdy(1) = -(1.0 - xi)  * inv2dy  ! SW
  dNdy(2) = -(1.0 + xi)  * inv2dy  ! SE
  dNdy(3) = +(1.0 - xi)  * inv2dy  ! NW
  dNdy(4) = +(1.0 + xi)  * inv2dy  ! NE

end subroutine smpm_grad

! ============================================================
!> sMPM bilinear shape functions with the same output interface as gimp_nodes.
!! Always returns ng=4 (SW, SE, NW, NE corners of the particle's home cell).
!! di_g and dj_g are in {0,1}: no extended support beyond the 4 corners.
subroutine smpm_nodes(xi, eta, dx, dy, Sw_g, dNdx_g, dNdy_g, di_g, dj_g, ng)
  real,    intent(in)  :: xi, eta       !< Particle local coords ∈ [-1,1] [nondim]
  real,    intent(in)  :: dx, dy        !< Cell dimensions [L]
  real,    intent(out) :: Sw_g(9)       !< Shape values (entries 1:4 set; 5:9 unused)
  real,    intent(out) :: dNdx_g(9)    !< ∂N/∂x [L-1]
  real,    intent(out) :: dNdy_g(9)    !< ∂N/∂y [L-1]
  integer, intent(out) :: di_g(9)      !< x-offsets from (ic-1): 0 or 1
  integer, intent(out) :: dj_g(9)      !< y-offsets from (jc-1): 0 or 1
  integer, intent(out) :: ng           !< Number of contributing nodes (always 4)
  real :: N(4), gx(4), gy(4)

  call smpm_shape(xi, eta, N)
  call smpm_grad(xi, eta, dx, dy, gx, gy)
  ng = 4
  ! Corner ordering matches smpm_shape: 1=SW, 2=SE, 3=NW, 4=NE
  Sw_g(1:4)   = N
  dNdx_g(1:4) = gx
  dNdy_g(1:4) = gy
  di_g(1:4) = [0, 1, 0, 1]
  dj_g(1:4) = [0, 0, 1, 1]
end subroutine smpm_nodes

! ============================================================
! Private helper: 1D GIMP shape function and gradient (Bardenhagen & Kober 2004).
!
! xpI: signed distance (particle minus node) [L > 0 means particle is right of node]
! lp:  particle half-width in this direction [L >= 0]
! h:   cell size [L > 0]
! Nval:   shape function value [nondim]
! dNval:  dN/d(x_particle) = dN/d(xpI) [L-1]
!
! When lp = 0 this reduces exactly to the standard linear hat function.
pure subroutine gimp1d(xpI, lp, h, Nval, dNval)
  real, intent(in)  :: xpI, lp, h
  real, intent(out) :: Nval, dNval
  real :: ax, tmp

  ax = abs(xpI)
  if (ax < lp) then
    ! Inner plateau: quadratic centred on node (only active when lp > 0)
    Nval  = 1.0 - (xpI*xpI + lp*lp) / (2.0 * h * lp)
    dNval = -xpI / (h * lp)
  elseif (ax < h - lp) then
    ! Linear ramp: same as standard sMPM
    Nval  = 1.0 - ax / h
    dNval = -sign(1.0, xpI) / h
  elseif (ax < h + lp) then
    ! Outer quadratic tail: extends support beyond standard cell boundary
    tmp   = h + lp - ax
    Nval  = tmp * tmp / (4.0 * h * lp)
    dNval = -sign(1.0, xpI) * tmp / (2.0 * h * lp)
  else
    Nval  = 0.0
    dNval = 0.0
  endif
end subroutine gimp1d

! ============================================================
!> GIMP (Generalized Interpolation MPM) shape functions and gradients.
!!
!! Computes N, dN/dx, dN/dy for all grid nodes whose GIMP support overlaps
!! particle p at local coordinate (xi, eta) ∈ [-1,1]² with physical half-widths
!! (lp_x, lp_y) [L] in a cell of size dx × dy [L].
!!
!! The standard 4-corner support (di,dj ∈ {0,1}²) is always included.
!! When lp_x > 0 or lp_y > 0, up to 5 extended nodes (di or dj ∈ {-1,2}) may
!! also contribute, giving at most 9 nodes total.
!!
!! Node index mapping for particle's home cell (ic, jc):
!!   I_Bu = ic - 1 + di_g(k),  di_g ∈ {-1, 0, 1, 2}
!!   J_Bu = jc - 1 + dj_g(k),  dj_g ∈ {-1, 0, 1, 2}
!!
!! Caller must bounds-check I_Bu and J_Bu against the allocated array extents.
subroutine gimp_nodes(xi, eta, lp_x, lp_y, dx, dy, Sw_g, dNdx_g, dNdy_g, di_g, dj_g, ng)
  real,    intent(in)  :: xi, eta       !< Particle local coords ∈ [-1,1]
  real,    intent(in)  :: lp_x, lp_y   !< Particle half-widths [L]
  real,    intent(in)  :: dx, dy       !< Cell dimensions [L]
  real,    intent(out) :: Sw_g(9)      !< GIMP shape values at non-zero nodes [nondim]
  real,    intent(out) :: dNdx_g(9)    !< ∂N/∂x [L-1]
  real,    intent(out) :: dNdy_g(9)    !< ∂N/∂y [L-1]
  integer, intent(out) :: di_g(9)      !< x-offsets from (ic-1), range -1..2
  integer, intent(out) :: dj_g(9)      !< y-offsets from (jc-1), range -1..2
  integer, intent(out) :: ng           !< Number of non-zero contributing nodes

  ! Full GIMP (Bardenhagen & Kober 2004 / Huth et al. 2021).
  ! 1D arrays for the 4 candidate x-positions (di = -1, 0, 1, 2 → kx = 0, 1, 2, 3)
  ! and 4 candidate y-positions (dj = -1, 0, 1, 2 → ky = 0, 1, 2, 3).
  ! Physical distance from node to particle for x-candidate kx:
  !   xpI = (xi - xi_I) * dx/2,  where xi_I = 2*di - 1 = 2*(kx-1) - 1
  !   so xpI = (xi - (2*(kx-1) - 1)) * dx/2 = (xi + 3 - 2*kx) * dx/2
  real :: Nx(0:3), dNx(0:3)
  real :: Ny(0:3), dNy(0:3)
  integer :: kx, ky, di_val, dj_val

  do kx = 0, 3
    call gimp1d((xi + 3.0 - 2.0*real(kx)) * 0.5 * dx, lp_x, dx, Nx(kx), dNx(kx))
  enddo
  do ky = 0, 3
    call gimp1d((eta + 3.0 - 2.0*real(ky)) * 0.5 * dy, lp_y, dy, Ny(ky), dNy(ky))
  enddo

  ng = 0
  do ky = 0, 3
    if (Ny(ky) == 0.0) cycle
    dj_val = ky - 1
    do kx = 0, 3
      if (Nx(kx) == 0.0) cycle
      di_val = kx - 1
      ng = ng + 1
      Sw_g(ng)    = Nx(kx) * Ny(ky)
      dNdx_g(ng) = dNx(kx) * Ny(ky)
      dNdy_g(ng) = Nx(kx) * dNy(ky)
      di_g(ng)   = di_val
      dj_g(ng)   = dj_val
    enddo
  enddo

end subroutine gimp_nodes

! ============================================================
! Private helper: swap all fields of particles pa and pb.
subroutine swap_particles(MPM, pa, pb)
  type(MPM_CS), intent(inout) :: MPM
  integer,      intent(in)    :: pa, pb

  real :: r ; integer :: k ; integer(kind=8) :: i8
  integer :: ii

#define SWAP_REAL(f)  r = MPM%f(pa); MPM%f(pa) = MPM%f(pb); MPM%f(pb) = r
#define SWAP_INT(f)   ii = MPM%f(pa); MPM%f(pa) = MPM%f(pb); MPM%f(pb) = ii
#define SWAP_INT8(f)  i8 = MPM%f(pa); MPM%f(pa) = MPM%f(pb); MPM%f(pb) = i8

  SWAP_REAL(xi) ; SWAP_REAL(eta)
  SWAP_INT(ci)  ; SWAP_INT(cj)
  SWAP_REAL(up) ; SWAP_REAL(vp)
  SWAP_REAL(up_g) ; SWAP_REAL(vp_g)
  SWAP_REAL(H) ; SWAP_REAL(PVolume) ; SWAP_REAL(GVolume)
  SWAP_REAL(Lx) ; SWAP_REAL(Ly) ; SWAP_REAL(Lx0) ; SWAP_REAL(Ly0)
  SWAP_REAL(strain_x) ; SWAP_REAL(strain_y)
  do k = 1, 4 ; r = MPM%Fdef(k,pa); MPM%Fdef(k,pa)=MPM%Fdef(k,pb); MPM%Fdef(k,pb)=r ; enddo
  do k = 1, 4 ; r = MPM%GradVel(k,pa); MPM%GradVel(k,pa)=MPM%GradVel(k,pb); MPM%GradVel(k,pb)=r ; enddo
  do k = 1, 2 ; r = MPM%GradZs(k,pa); MPM%GradZs(k,pa)=MPM%GradZs(k,pb); MPM%GradZs(k,pb)=r ; enddo
  do k = 1, 2 ; r = MPM%GradH(k,pa);  MPM%GradH(k,pa) =MPM%GradH(k,pb);  MPM%GradH(k,pb) =r ; enddo
  SWAP_REAL(AGLen) ; SWAP_REAL(eta_visc) ; SWAP_REAL(smb)
  SWAP_INT(status) ; SWAP_INT8(pid)

#undef SWAP_REAL
#undef SWAP_INT
#undef SWAP_INT8

end subroutine swap_particles

! ============================================================
! Private helper: in-place resize of a real rank-1 allocatable.
subroutine grow_real_1d(arr, n_new)
  real, allocatable, intent(inout) :: arr(:)
  integer, intent(in) :: n_new
  real, allocatable :: tmp(:)
  integer :: n_old
  n_old = size(arr)
  allocate(tmp(n_new), source=0.0)
  tmp(1:n_old) = arr(1:n_old)
  call move_alloc(tmp, arr)
end subroutine grow_real_1d

subroutine grow_real_2d(arr, m, n_new)
  real, allocatable, intent(inout) :: arr(:,:)
  integer, intent(in) :: m, n_new
  real, allocatable :: tmp(:,:)
  integer :: n_old
  n_old = size(arr, 2)
  allocate(tmp(m, n_new), source=0.0)
  tmp(:, 1:n_old) = arr(:, 1:n_old)
  call move_alloc(tmp, arr)
end subroutine grow_real_2d

subroutine grow_int_1d(arr, n_new)
  integer, allocatable, intent(inout) :: arr(:)
  integer, intent(in) :: n_new
  integer, allocatable :: tmp(:)
  integer :: n_old
  n_old = size(arr)
  allocate(tmp(n_new), source=0)
  tmp(1:n_old) = arr(1:n_old)
  call move_alloc(tmp, arr)
end subroutine grow_int_1d

subroutine grow_int8_1d(arr, n_new)
  integer(kind=8), allocatable, intent(inout) :: arr(:)
  integer, intent(in) :: n_new
  integer(kind=8), allocatable :: tmp(:)
  integer :: n_old
  n_old = size(arr)
  allocate(tmp(n_new), source=0_8)
  tmp(1:n_old) = arr(1:n_old)
  call move_alloc(tmp, arr)
end subroutine grow_int8_1d

! ============================================================
! Private helper: integer to string (for log messages).
function int_to_str(n) result(s)
  integer, intent(in) :: n
  character(len=20) :: s
  write(s,'(i20)') n
  s = adjustl(s)
end function int_to_str

end module MOM_ice_shelf_MPM
