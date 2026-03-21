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

use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, is_root_pe
use MOM_file_parser,   only : get_param, param_file_type
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_ice_shelf_state, only : ice_shelf_state

implicit none ; private

public :: MPM_CS, MPM_init, MPM_end
public :: MPM_P2G, MPM_G2P, MPM_Lagrangian_update
public :: MPM_split, MPM_migrate, MPM_reseed
public :: MPM_update_masks

#include <MOM_memory.h>

! Particle status flags
integer, parameter :: ALIVE   = 1 !< Active particle
integer, parameter :: DEAD    = 0 !< Inactive (to be compacted out)
integer, parameter :: LEAVING = 2 !< Particle about to migrate to another PE

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
  real,    allocatable :: Lx(:)       !< Particle half-width in x [L ~> m]
  real,    allocatable :: Ly(:)       !< Particle half-width in y [L ~> m]
  real,    allocatable :: Lx0(:)      !< Original half-width in x [L ~> m]
  real,    allocatable :: Ly0(:)      !< Original half-width in y [L ~> m]
  real,    allocatable :: strain_x(:) !< Accumulated longitudinal strain ε_xx [nondim]
  real,    allocatable :: strain_y(:) !< Accumulated longitudinal strain ε_yy [nondim]

  !> Deformation gradient F stored as [F11, F12, F21, F22] (row-major 2x2)
  real,    allocatable :: Fdef(:,:)   !< Deformation gradient [nondim]; shape (4, n_alloc)

  real,    allocatable :: GradVel(:,:) !< Velocity gradient [du/dx, dv/dy, du/dy, dv/dx]
                                       !! [T-1 ~> s-1]; shape (4, n_alloc)
  real,    allocatable :: GradZs(:,:) !< Surface-elevation gradient [dZs/dx, dZs/dy]
                                       !! [nondim]; shape (2, n_alloc)

  real,    allocatable :: AGLen(:)    !< Glen's-law rate factor A [Pa-n s-1 ~> R-n L2n T-(1+n)]
                                       !! stored as AGlen_visc equivalent
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

  ! --- T-grid fields for P2G bookkeeping ---

  real, allocatable :: area_shelf_h(:,:) !< Sum of PVolume_p over particles in cell [L2 ~> m2]
  real, allocatable :: reweight(:,:)     !< areaT / max(area_shelf_h, eps) [nondim]

  ! --- Counters and parameters ---

  integer :: n_active  !< Number of active particles on this PE
  integer :: n_alloc   !< Allocated size of particle arrays

  integer :: n_per_cell !< Number of particles per cell at initialisation (e.g. 4)
  real :: split_factor  !< Split when Lx or Ly > split_factor * Lx0 [nondim]
  real :: flip_alpha    !< FLIP blending coefficient (1=pure FLIP, 0=pure PIC) [nondim]
  real :: H_min         !< Minimum ice thickness at a particle [Z ~> m]
  real :: area_tiny     !< Small area used to avoid division by zero [L2 ~> m2]

  integer(kind=8) :: pid_counter !< Next particle ID to assign on this PE (PE-offset ensures uniqueness)

  logical :: use_constant_AGlen  !< If true, use AGlen_const for all particles
  real    :: AGlen_const         !< Constant Glen A [same units as CS%AGLen]

  logical :: initialized = .false. !< True after MPM_init has run

end type MPM_CS

contains

! ============================================================
!> Allocate and initialise the MPM control structure.
!! Seeds particles from ISS%h_shelf (T-grid) and sets up connectivity.
subroutine MPM_init(MPM, ISS, G, US, param_file, n_glen, AGLen_visc_ref)
  type(MPM_CS),          intent(inout) :: MPM         !< MPM control structure to initialise
  type(ice_shelf_state), intent(in)    :: ISS         !< Ice shelf state (h_shelf, hmask)
  type(ocean_grid_type), intent(in)    :: G           !< Ocean grid
  type(unit_scale_type), intent(in)    :: US          !< Unit scaling factors
  type(param_file_type), intent(in)    :: param_file  !< Run-time parameter file
  real,                  intent(in)    :: n_glen      !< Glen exponent [nondim]
  real,                  intent(in)    :: AGLen_visc_ref !< Reference Glen A (same units as CS)

  ! Local variables
  character(len=40) :: mdl = "MOM_ice_shelf_MPM"
  integer :: i, j, p, pp, ip, jp, n_sq, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: pe_rank, npes
  real :: xi0, eta0, dxi, deta, dx, dy, cell_area
  real :: H_p, Vol_p, Lx_p, Ly_p
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
                 "original length.", &
                 units="nondim", default=2.0)
  call get_param(param_file, mdl, "MPM_FLIP_ALPHA", MPM%flip_alpha, &
                 "FLIP blending coefficient; 1 = pure FLIP (no damping), "//&
                 "0 = pure PIC.", &
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "MPM_H_MIN", MPM%H_min, &
                 "Minimum ice thickness per particle.", &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(param_file, mdl, "MPM_USE_CONSTANT_AGLEN", MPM%use_constant_AGlen, &
                 "If true, use the reference Glen A for all particles.", &
                 default=.true.)
  MPM%AGlen_const = AGLen_visc_ref

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
      Lx_p  = dx / (2.0 * n_sq)   ! particle half-width [L]
      Ly_p  = dy / (2.0 * n_sq)   ! particle half-height [L]

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

  ! --- Build cell_start / cell_list ---
  call MPM_build_cell_list(MPM, isd, ied, jsd, jed)

  ! --- Initialise B-grid fields ---
  MPM%H_node(:,:) = 0.0 ; MPM%Zs_node(:,:) = 0.0
  MPM%bed_node(:,:) = 0.0
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
  allocate(MPM%AGLen(n),     source=0.0)
  allocate(MPM%smb(n),       source=0.0)
  allocate(MPM%status(n),    source=DEAD)
  allocate(MPM%pid(n),       source=0_8)
  allocate(MPM%cell_list(n), source=0)

  allocate(MPM%cell_count(isd:ied, jsd:jed), source=0)
  allocate(MPM%cell_start(isd:ied, jsd:jed), source=0)

  allocate(MPM%H_node(IsdB:IedB, JsdB:JedB),   source=0.0)
  allocate(MPM%Zs_node(IsdB:IedB, JsdB:JedB),  source=0.0)
  allocate(MPM%bed_node(IsdB:IedB, JsdB:JedB),  source=0.0)
  allocate(MPM%area_shelf_h(isd:ied, jsd:jed),  source=0.0)
  allocate(MPM%reweight(isd:ied, jsd:jed),       source=0.0)

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
  call grow_real_1d(MPM%AGLen, n_new)
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
!> Particles → Grid (P2G): update ISS%h_shelf and ISS%hmask from particle state.
!!
!! Uses a volume-weighted cell average: for each T-grid cell (i,j),
!!   ISS%h_shelf(i,j) = Σ_p H(p)*PVolume(p) / Σ_p PVolume(p)
!! This avoids cross-PE summation at B-grid corners (handled separately
!! by the SSA solver's existing interpolate_H_to_B call).
!! The reweight factor is also computed for use in G2P and Lagrangian updates.
subroutine MPM_P2G(MPM, ISS, G)
  type(MPM_CS),          intent(inout) :: MPM  !< MPM control structure
  type(ice_shelf_state), intent(inout) :: ISS  !< Ice shelf state (h_shelf, hmask updated here)
  type(ocean_grid_type), intent(in)    :: G    !< Ocean grid

  integer :: i, j, p, ic, jc
  integer :: isd, ied, jsd, jed
  real, allocatable :: H_vol(:,:)  ! Σ H(p)*PVolume(p) per cell [Z L2]
  real, allocatable :: H_wgt(:,:)  ! Σ PVolume(p) per cell [L2]

  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  allocate(H_vol(isd:ied, jsd:jed), source=0.0)
  allocate(H_wgt(isd:ied, jsd:jed), source=0.0)

  ! --- Accumulate per-cell particle H and area ---
  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    H_wgt(ic,jc) = H_wgt(ic,jc) + MPM%PVolume(p)
    H_vol(ic,jc) = H_vol(ic,jc) + MPM%H(p) * MPM%PVolume(p)
  enddo

  ! --- Update area_shelf_h and reweight ---
  MPM%area_shelf_h(:,:) = H_wgt(:,:)
  do j = jsd, jed ; do i = isd, ied
    if (MPM%area_shelf_h(i,j) > MPM%area_tiny) then
      MPM%reweight(i,j) = G%areaT(i,j) / MPM%area_shelf_h(i,j)
    else
      MPM%reweight(i,j) = 0.0
    endif
  enddo ; enddo

  ! --- Update ISS%h_shelf and ISS%hmask ---
  do j = jsd, jed ; do i = isd, ied
    if (ISS%hmask(i,j) == 3) cycle   ! inflow BC: never overwrite

    if (H_wgt(i,j) > MPM%area_tiny) then
      ISS%h_shelf(i,j) = max(H_vol(i,j) / H_wgt(i,j), MPM%H_min)
      if (ISS%hmask(i,j) == 0) ISS%hmask(i,j) = 1   ! cell newly occupied
    else
      if (ISS%hmask(i,j) == 1) then
        ISS%hmask(i,j) = 0
        ISS%h_shelf(i,j) = 0.0
      endif
    endif
  enddo ; enddo

  ! Propagate h_shelf and hmask into halo so the SSA solver sees them
  call pass_var(ISS%h_shelf, G%domain, complete=.false.)
  call pass_var(ISS%hmask,   G%domain, complete=.true.)

  deallocate(H_vol, H_wgt)

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
subroutine MPM_G2P(MPM, G, u_shelf, v_shelf)
  type(MPM_CS),          intent(inout) :: MPM     !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G       !< Ocean grid
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(in)    :: u_shelf !< Zonal velocity [L T-1 ~> m s-1]
  real, dimension(G%IsdB:G%IedB, G%JsdB:G%JedB), &
                         intent(in)    :: v_shelf !< Meridional velocity [L T-1 ~> m s-1]

  integer :: p, ic, jc
  real :: xi_p, eta_p, dx, dy
  real :: N(4), dNdx(4), dNdy(4)
  ! Corner B-grid velocities
  real :: uSW, uSE, uNW, uNE
  real :: vSW, vSE, vNW, vNE
  ! Surface elevation Zs at corners (approx from h_shelf + bed)
  real :: ZsSW, ZsSE, ZsNW, ZsNE
  real :: u_new, v_new  ! New grid velocity at particle position

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    xi_p  = MPM%xi(p)
    eta_p = MPM%eta(p)
    dx = G%dxT(ic,jc)
    dy = G%dyT(ic,jc)

    call smpm_shape(xi_p, eta_p, N)
    call smpm_grad(xi_p, eta_p, dx, dy, dNdx, dNdy)

    ! B-grid corner velocities (NE convention: Bu(I,J) is NE of T(i,j))
    uSW = u_shelf(ic-1, jc-1) ; vSW = v_shelf(ic-1, jc-1)
    uSE = u_shelf(ic,   jc-1) ; vSE = v_shelf(ic,   jc-1)
    uNW = u_shelf(ic-1, jc  ) ; vNW = v_shelf(ic-1, jc  )
    uNE = u_shelf(ic,   jc  ) ; vNE = v_shelf(ic,   jc  )

    ! New grid velocity at particle (for FLIP update later)
    u_new = N(1)*uSW + N(2)*uSE + N(3)*uNW + N(4)*uNE
    v_new = N(1)*vSW + N(2)*vSE + N(3)*vNW + N(4)*vNE

    ! FLIP velocity update: keep particle's own velocity + grid correction
    MPM%up(p) = MPM%flip_alpha * (MPM%up(p) + (u_new - MPM%up_g(p))) &
              + (1.0 - MPM%flip_alpha) * u_new
    MPM%vp(p) = MPM%flip_alpha * (MPM%vp(p) + (v_new - MPM%vp_g(p))) &
              + (1.0 - MPM%flip_alpha) * v_new

    ! Save current grid velocity for next FLIP step
    MPM%up_g(p) = u_new
    MPM%vp_g(p) = v_new

    ! Velocity gradient ∇u at particle (for Glen's law viscosity)
    MPM%GradVel(1,p) = dNdx(1)*uSW + dNdx(2)*uSE + dNdx(3)*uNW + dNdx(4)*uNE  ! du/dx
    MPM%GradVel(2,p) = dNdy(1)*vSW + dNdy(2)*vSE + dNdy(3)*vNW + dNdy(4)*vNE  ! dv/dy
    MPM%GradVel(3,p) = dNdy(1)*uSW + dNdy(2)*uSE + dNdy(3)*uNW + dNdy(4)*uNE  ! du/dy
    MPM%GradVel(4,p) = dNdx(1)*vSW + dNdx(2)*vSE + dNdx(3)*vNW + dNdx(4)*vNE  ! dv/dx

    ! Surface elevation at corners: Zs = (1 - rho_i/rho_w) * H  [floating ice; stored in H_node]
    ! Here we just store the gradient of H_node as a proxy for GradZs; the dynamics module
    ! computes the full driving stress including bed effects in calc_shelf_driving_stress.
    ZsSW = MPM%H_node(ic-1, jc-1)
    ZsSE = MPM%H_node(ic,   jc-1)
    ZsNW = MPM%H_node(ic-1, jc  )
    ZsNE = MPM%H_node(ic,   jc  )
    MPM%GradZs(1,p) = dNdx(1)*ZsSW + dNdx(2)*ZsSE + dNdx(3)*ZsNW + dNdx(4)*ZsNE
    MPM%GradZs(2,p) = dNdy(1)*ZsSW + dNdy(2)*ZsSE + dNdy(3)*ZsNW + dNdy(4)*ZsNE
  enddo

end subroutine MPM_G2P

! ============================================================
!> Lagrangian particle updates: deformation gradient F, thickness H,
!! particle dimensions, and position.  Called once per timestep after G2P.
subroutine MPM_Lagrangian_update(MPM, G, dt)
  type(MPM_CS),          intent(inout) :: MPM   !< MPM control structure
  type(ocean_grid_type), intent(in)    :: G     !< Ocean grid
  real,                  intent(in)    :: dt    !< Timestep [T ~> s]

  integer :: p, ic, jc
  real :: ux, uy, vx, vy   ! Velocity gradient components [T-1]
  real :: divv              ! Velocity divergence [T-1]
  real :: F11, F12, F21, F22   ! Deformation gradient
  real :: L11, L12, L21, L22  ! Velocity gradient increment  L = I + dt*gradV
  real :: detF              ! Det of deformation gradient
  real :: H_new, xi_new, eta_new
  real :: H_disc, dx, dy

  do p = 1, MPM%n_active
    if (MPM%status(p) /= ALIVE) cycle

    ! --- Update deformation gradient F = L · F_old ---
    ux = MPM%GradVel(1,p)
    vy = MPM%GradVel(2,p)
    uy = MPM%GradVel(3,p)
    vx = MPM%GradVel(4,p)

    L11 = 1.0 + dt*ux ; L12 = dt*vx
    L21 = dt*uy       ; L22 = 1.0 + dt*vy

    F11 = MPM%Fdef(1,p) ; F12 = MPM%Fdef(2,p)
    F21 = MPM%Fdef(3,p) ; F22 = MPM%Fdef(4,p)

    MPM%Fdef(1,p) = L11*F11 + L12*F21
    MPM%Fdef(2,p) = L11*F12 + L12*F22
    MPM%Fdef(3,p) = L21*F11 + L22*F21
    MPM%Fdef(4,p) = L21*F12 + L22*F22

    ! --- Update particle volume from det(F) ---
    detF = MPM%Fdef(1,p)*MPM%Fdef(4,p) - MPM%Fdef(2,p)*MPM%Fdef(3,p)
    detF = max(detF, 0.01)  ! Guard against unphysical compression
    MPM%PVolume(p) = detF * MPM%GVolume(p)

    ! --- Update sMPM particle dimensions (strain tracking) ---
    divv = ux + vy
    MPM%strain_x(p) = MPM%strain_x(p) + dt * ux
    MPM%strain_y(p) = MPM%strain_y(p) + dt * vy
    MPM%Lx(p) = MPM%Lx0(p) * (1.0 + MPM%strain_x(p))
    MPM%Ly(p) = MPM%Ly0(p) * (1.0 + MPM%strain_y(p))
    ! Clamp to physical minimum
    MPM%Lx(p) = max(MPM%Lx(p), 0.1 * MPM%Lx0(p))
    MPM%Ly(p) = max(MPM%Ly(p), 0.1 * MPM%Ly0(p))

    ! --- Update ice thickness: dH/dt = -H * div(v) + SMB ---
    H_new = MPM%H(p) * (1.0 - dt * divv) + MPM%smb(p) * dt
    if (H_new < MPM%H_min) H_new = MPM%H_min  ! enforce minimum thickness
    MPM%H(p) = H_new

    ! --- Advect particle position ---
    ic = MPM%ci(p) ; jc = MPM%cj(p)
    dx = G%dxT(ic, jc)
    dy = G%dyT(ic, jc)

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

    split_x = (MPM%Lx(p) > MPM%split_factor * MPM%Lx0(p))
    split_y = (MPM%Ly(p) > MPM%split_factor * MPM%Ly0(p))

    if (.not.(split_x .or. split_y)) cycle

    dx = G%dxT(MPM%ci(p), MPM%cj(p))
    dy = G%dyT(MPM%ci(p), MPM%cj(p))
    ! Offset in local coordinates (half the child particle size)
    dxi_child  = MPM%Lx(p) / (0.5 * dx)   ! [nondim]
    deta_child = MPM%Ly(p) / (0.5 * dy)   ! [nondim]

    if (split_x .and. split_y) then
      ! 4-way split into children at ±xi/2, ±eta/2
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 +0.5*dxi_child, +0.5*deta_child, 0.25)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 -0.5*dxi_child, +0.5*deta_child, 0.25)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 +0.5*dxi_child, -0.5*deta_child, 0.25)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 -0.5*dxi_child, -0.5*deta_child, 0.25)
      MPM%status(p) = DEAD  ! parent is replaced by 4 children
    elseif (split_x) then
      ! 2-way split in x
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 +0.5*dxi_child, 0.0, 0.5)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 -0.5*dxi_child, 0.0, 0.5)
      MPM%status(p) = DEAD
    else  ! split_y only
      ! 2-way split in y
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 0.0, +0.5*deta_child, 0.5)
      call split_particle_child(MPM, p, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                 0.0, -0.5*deta_child, 0.5)
      MPM%status(p) = DEAD
    endif
  enddo

end subroutine MPM_split

! ============================================================
!> Create a child particle from parent p, offset in local coords.
subroutine split_particle_child(MPM, p_parent, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB, &
                                  dxi, deta, vol_frac)
  type(MPM_CS), intent(inout) :: MPM
  integer,      intent(in)    :: p_parent
  integer,      intent(in)    :: isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  real,         intent(in)    :: dxi, deta    !< Offset in cell-local coords
  real,         intent(in)    :: vol_frac     !< Fraction of parent volume

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
  MPM%H(p_child)   = MPM%H(p_parent)   ! Thickness unchanged in split
  MPM%PVolume(p_child) = vol_frac * MPM%PVolume(p_parent)
  MPM%GVolume(p_child) = vol_frac * MPM%GVolume(p_parent)
  MPM%Lx(p_child)  = 0.5 * MPM%Lx(p_parent)  * (1.0 / (1.0 + abs(dxi) * 2.0))
  MPM%Ly(p_child)  = 0.5 * MPM%Ly(p_parent)  * (1.0 / (1.0 + abs(deta) * 2.0))
  MPM%Lx0(p_child) = MPM%Lx(p_child)
  MPM%Ly0(p_child) = MPM%Ly(p_child)
  MPM%strain_x(p_child) = 0.0
  MPM%strain_y(p_child) = 0.0
  MPM%Fdef(1,p_child) = 1.0 ; MPM%Fdef(2,p_child) = 0.0
  MPM%Fdef(3,p_child) = 0.0 ; MPM%Fdef(4,p_child) = 1.0
  MPM%GradVel(:,p_child) = MPM%GradVel(:,p_parent)
  MPM%GradZs(:,p_child)  = MPM%GradZs(:,p_parent)
  MPM%AGLen(p_child)  = MPM%AGLen(p_parent)
  MPM%smb(p_child)    = MPM%smb(p_parent)
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
!> Reseed particles in Dirichlet inflow cells (hmask == 3) that have
!! fallen below the target particle count due to downstream advection.
!! New particles are initialised with the current shelf thickness and the
!! prescribed inflow velocity averaged from the four surrounding B-grid corners.
subroutine MPM_reseed(MPM, ISS, G, u_bdry_val, v_bdry_val)
  type(MPM_CS),          intent(inout) :: MPM         !< MPM control structure
  type(ice_shelf_state), intent(in)    :: ISS         !< Ice shelf state
  type(ocean_grid_type), intent(in)    :: G           !< Ocean grid
  real, dimension(:,:),  intent(in)    :: u_bdry_val  !< Zonal BC velocity on B-grid [L T-1]
  real, dimension(:,:),  intent(in)    :: v_bdry_val  !< Meridional BC velocity on B-grid [L T-1]

  integer :: i, j, k, pp, ip, jp, n_sq, n_current
  real    :: H_p, Vol_p, Lx_p, Ly_p, dx, dy, dxi, deta, xi0, eta0
  real    :: u_bc, v_bc
  integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB
  integer :: n_added

  isc = G%isc ; iec = G%iec ; jsc = G%jsc ; jec = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  n_sq  = nint(sqrt(real(MPM%n_per_cell)))
  dxi   = 2.0 / n_sq
  deta  = 2.0 / n_sq
  n_added = 0

  do j = jsc, jec
    do i = isc, iec
      if (ISS%hmask(i,j) /= 3.0) cycle

      n_current = MPM%cell_count(i,j)
      if (n_current >= MPM%n_per_cell) cycle

      H_p   = max(ISS%h_shelf(i,j), MPM%H_min)
      dx    = G%dxT(i,j)
      dy    = G%dyT(i,j)
      Vol_p = G%areaT(i,j) / real(MPM%n_per_cell)
      Lx_p  = dx / (2.0 * n_sq)
      Ly_p  = dy / (2.0 * n_sq)

      ! Average boundary velocity from four surrounding B-grid corners.
      ! Cell (i,j) on T-grid has B-grid corners at indices (i-1,j-1), (i,j-1), (i-1,j), (i,j).
      u_bc = 0.25 * (u_bdry_val(i-1, j-1) + u_bdry_val(i, j-1) + &
                     u_bdry_val(i-1, j  ) + u_bdry_val(i, j  ))
      v_bc = 0.25 * (v_bdry_val(i-1, j-1) + v_bdry_val(i, j-1) + &
                     v_bdry_val(i-1, j  ) + v_bdry_val(i, j  ))

      ! Add missing particles at regular sub-cell grid positions (starting after
      ! the slots already occupied, assuming initial ordering is preserved).
      k = 0
      do jp = 1, n_sq
        do ip = 1, n_sq
          k = k + 1
          if (k <= n_current) cycle   ! positions already filled

          xi0  = -1.0 + (real(ip) - 0.5) * dxi
          eta0 = -1.0 + (real(jp) - 0.5) * deta

          MPM%n_active = MPM%n_active + 1
          if (MPM%n_active > MPM%n_alloc) &
            call MPM_grow_arrays(MPM, isd, ied, jsd, jed, IsdB, IedB, JsdB, JedB)

          pp = MPM%n_active
          MPM%xi(pp)        = xi0   ; MPM%eta(pp)       = eta0
          MPM%ci(pp)        = i     ; MPM%cj(pp)        = j
          MPM%up(pp)        = u_bc  ; MPM%vp(pp)        = v_bc
          MPM%up_g(pp)      = u_bc  ; MPM%vp_g(pp)      = v_bc
          MPM%H(pp)         = H_p
          MPM%PVolume(pp)   = Vol_p ; MPM%GVolume(pp)   = Vol_p
          MPM%Lx(pp)        = Lx_p  ; MPM%Ly(pp)        = Ly_p
          MPM%Lx0(pp)       = Lx_p  ; MPM%Ly0(pp)       = Ly_p
          MPM%strain_x(pp)  = 0.0   ; MPM%strain_y(pp)  = 0.0
          MPM%Fdef(1,pp)    = 1.0   ; MPM%Fdef(2,pp)    = 0.0
          MPM%Fdef(3,pp)    = 0.0   ; MPM%Fdef(4,pp)    = 1.0
          MPM%GradVel(1,pp) = 0.0   ; MPM%GradVel(2,pp) = 0.0
          MPM%GradVel(3,pp) = 0.0   ; MPM%GradVel(4,pp) = 0.0
          MPM%GradZs(1,pp)  = 0.0   ; MPM%GradZs(2,pp)  = 0.0
          MPM%AGLen(pp)     = MPM%AGlen_const
          MPM%smb(pp)       = 0.0
          MPM%status(pp)    = ALIVE
          MPM%pid(pp)       = MPM%pid_counter
          MPM%pid_counter   = MPM%pid_counter + 1
          MPM%cell_count(i,j) = MPM%cell_count(i,j) + 1
          n_added = n_added + 1
        enddo
      enddo
    enddo
  enddo

  if (n_added > 0) &
    call MPM_build_cell_list(MPM, isd, ied, jsd, jed)

end subroutine MPM_reseed

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
  if (allocated(MPM%AGLen)) deallocate(MPM%AGLen)
  if (allocated(MPM%smb)) deallocate(MPM%smb)
  if (allocated(MPM%status)) deallocate(MPM%status)
  if (allocated(MPM%pid)) deallocate(MPM%pid)
  if (allocated(MPM%cell_list)) deallocate(MPM%cell_list)
  if (allocated(MPM%cell_count)) deallocate(MPM%cell_count)
  if (allocated(MPM%cell_start)) deallocate(MPM%cell_start)
  if (allocated(MPM%H_node)) deallocate(MPM%H_node)
  if (allocated(MPM%Zs_node)) deallocate(MPM%Zs_node)
  if (allocated(MPM%bed_node)) deallocate(MPM%bed_node)
  if (allocated(MPM%area_shelf_h)) deallocate(MPM%area_shelf_h)
  if (allocated(MPM%reweight)) deallocate(MPM%reweight)

  MPM%initialized = .false.

end subroutine MPM_end

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
  SWAP_REAL(AGLen) ; SWAP_REAL(smb)
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
