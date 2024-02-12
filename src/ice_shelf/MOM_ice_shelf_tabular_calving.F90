!> Tabular calving routines for MOM6-IS.
module MOM_ice_shelf_tabular_calving

  ! This file is part of MOM6. See LICENSE.md for the license.

use MOM_array_transform, only : rotate_array
use mpp_mod, only : mpp_npes, mpp_pe, mpp_root_pe, NULL_PE
use MOM_coms, only : max_across_PEs, min_across_PEs
!use mpp_domains_mod, only: mpp_get_neighbor_pe, NORTH, SOUTH, EAST, WEST
use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock, only : CLOCK_COMPONENT, CLOCK_ROUTINE
use MOM_dyn_horgrid, only : dyn_horgrid_type, create_dyn_horgrid, destroy_dyn_horgrid
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_grid, only : MOM_grid_init, ocean_grid_type
use MOM_get_input, only : directories, Get_MOM_input
use MOM_coms, only : reproducing_sum
use MOM_checksums, only : hchksum, qchksum, chksum, uchksum, vchksum, uvchksum
use mpp_mod, only : mpp_sync_self
use MOM_file_parser, only : read_param, get_param, log_param, log_version, param_file_type
use MOM_io, only : slasher, file_exists
use MOM_interpolate, only : init_external_field, time_interp_external, time_interp_external_init
use MOM_interpolate, only : external_field
use MOM_unit_scaling, only : unit_scale_type
use MOM_ice_shelf_state, only : ice_shelf_state
use MOM_time_manager, only : time_type

  implicit none ; private

public initialize_tabular_calving, tabular_calving_end

!> Structure that describes the ice shelf calving state
type, public :: tabular_calving_state
  real, pointer, dimension(:,:) :: &
    tabular_calve_mask => NULL() !< Mask used to indicate cells ready to be calved and
                                 !! converted to bonded-particle tabular icebergs
                                 !! 0: not ready to calve
                                 !! >1: ready to calve

  logical :: tabular_calving_from_file !< If true, read in tabular calving mask from a file

  logical :: debug             !< If true, write verbose output

  type(external_field) :: calving_mask_handle !< Handle for reading the time-interpolated
                                              !! ice shelf tabular calving mask from a file
end type tabular_calving_state

contains

!> Initializes ice shelf tabular calving data, parameters, and diagnostics
subroutine initialize_tabular_calving(param_file, TC, G)
  type(param_file_type),        intent(in)    :: param_file !< A structure to parse for run-time parameters
  type(tabular_calving_state), pointer    :: TC !< A pointer to the tabular calving structure
  type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.
  character(len=40) :: mdl = "MOM_ice_shelf_tabular_calving" ! This module's name
  character(len=240) :: inputdir, TC_file, filename
  character(len=120) :: TC_mask_var  ! The name of shelf mass in the file.
  integer :: isd, ied, jsd, jed
  isd = G%isd ; jsd = G%jsd ; ied = G%ied ; jed = G%jed

  if (associated(TC)) then
    call MOM_error(FATAL, "MOM_ice_shelf_calving.F90, initialize_tabular_calving: "// &
                          "called with an associated tabular_calving pointer.")
    return
  endif

  allocate(TC)
  !allocate(TC%c_id(isd:ied,jsd:jed), source=0.0 )
  allocate(TC%tabular_calve_mask(isd:ied,jsd:jed), source=0.0 )

  call get_param(param_file, mdl, "DEBUG", TC%debug, default=.false.)
  !call get_param(param_file, mdl, "TABULAR_BERG_PARTICLE_RADIUS", CS%tabular_rad, &
  !               "particle radius for iKID particles that calve from ice shelves",&
  !               units="m", default=1000.0, scale=US%m_to_L)

  !In lieu of a tabular calving law or rifting model, we can currently only read in the tabular
  !calving mask from a file
  call get_param(param_file, mdl, "TABULAR_CALVING_MASK_FROM_FILE", &
                 TC%tabular_calving_from_file, "Read the mass of the "//&
                 "ice shelf (every time step) from a file.", default=.true.)

  if (TC%tabular_calving_from_file) then

    call time_interp_external_init()

    call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")

    inputdir = slasher(inputdir)

    call get_param(param_file, mdl, "TABULAR_CALVING_MASK_FILE", TC_file, &
            "If TABULAR_CALVING_MASK_FROM_FILE = True, this is the file from "//&
            "which to read the tabular calving mask.", &
            default="tabular_calving_mask.nc")

    call get_param(param_file, mdl, "TABULAR_CALVING_MASK_VAR", TC_mask_var, &
               "The variable in TABULAR_CALVING_MASK_FILE with the tabular calving mask.", &
               default="tabular_calving_mask")

    filename = trim(slasher(inputdir))//trim(TC_file)
    call log_param(param_file, mdl, "INPUTDIR/TABULAR_CALVING_MASK_FILE", filename)

    if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
         "initialize_tabular_calving: Unable to open "//trim(filename))

    TC%calving_mask_handle = init_external_field(filename, TC_mask_var, &
                                                 MOM_domain=G%Domain, verbose=TC%debug)
  endif

end subroutine initialize_tabular_calving

!> Deallocates all memory associated with this module
subroutine tabular_calving_end(TC)
  type(tabular_calving_state), pointer :: TC !< A pointer to the ice shelf calving structure

  if (.not.associated(TC)) return

!  if (allocated(TC%berg_list)) deallocate(TC%berg_list)
  deallocate(TC%tabular_calve_mask)
!  deallocate(TC%c_id)
  deallocate(TC)
end subroutine tabular_calving_end

!!$!> Initialize iKID icebergs from a tabular calving mask.
!!$subroutine process_tabular_calving(G, CS, ISS, TC, frac_cberg_calved, Time)
!!$  ! Arguments
!!$  type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.
!!$  type(ice_shelf_CS),    pointer    :: CS   !< A pointer to the control structure returned
!!$                                            !! by a previous call to initialize_ice_shelf.
!!$  type(ice_shelf_state), pointer :: ISS     !< A structure with elements that describe
!!$                                            !! the ice-shelf state
!!$  type(tabular_calving_state), pointer :: TC !< A pointer to the tabular calving structure
!!$  real, dimension(:,:), intent(inout) :: frac_cberg_calved !< cell fraction of fully-calved bonded bergs
!!$                                                           !! from the ice sheet [nondim]
!!$  type(time_type),       intent(in)    :: Time !< The current model time
!!$  ! Local variables
!!$  integer :: max_TC_mask
!!$  integer :: i, j, is, ie, js, je, m, n
!!$
!!$  !TC%c_id(:,:) = 0
!!$
!!$  !First, remove any ice that has fully calved
!!$  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
!!$  do j=js,je ; do i=is,ie
!!$    if (frac_cberg_calved(i,j)>=1) then
!!$      TC%tabular_calve_mask(i,j)    = 0
!!$      ISS%mass_shelf(i,j)   = 0
!!$      ISS%area_shelf_h(i,j) = 0
!!$      ISS%hmask(i,j)        = 0
!!$      ISS%h_shelf(i,j)      = 0
!!$    elseif (frac_cberg_calved(i,j)>0 .and. frac_cberg_calved(i,j)<1) then
!!$      ISS%mass_shelf(i,j)   = (1-frac_cberg_calved(i,j))*ISS%mass_shelf(i,j)
!!$      ISS%area_shelf_h(i,j) = (1-frac_cberg_calved(i,j))*ISS%area_shelf_h(i,j)
!!$      ISS%hmask(i,j)=2
!!$    endif
!!$    !reset frac_cberg_calved
!!$    frac_cberg_calved(i,j)=0
!!$  enddo; enddo
!!$
!!$  call pass_var(ISS%mass_shelf,   G%domain, complete=.false.)
!!$  call pass_var(ISS%area_shelf_h, G%domain, complete=.false.)
!!$  call pass_var(ISS%hmask,        G%domain, complete=.false.)
!!$  call pass_var(ISS%h_shelf,      G%domain, complete=.true.)
!!$
!!$  !First, update the tabular calving mask
!!$  call update_tabular_calving_mask(G, CS, TC, Time)
!!$
!!$  ! max_TC_mask = maxval(TC%tabular_calve_mask)
!!$  ! call max_across_PEs(max_TC_mask)
!!$  ! if (max_TC_mask<=0) then !no calving
!!$  !   if (allocated(TC%berg_list)) deallocate(TC%berg_list)
!!$  !   return
!!$  ! endif
!!$
!!$  ! !We want to initialize (calve)  bonded-particle tabular bergs that cover the "berg cells" defined
!!$  ! !where mask>0. There may be multiple tabular bergs to calve in the domain,
!!$  ! !which may overlap multiple cells and PEs. The approach is to create a field of unique IDs
!!$  ! !that are each associated with a different tabular berg. This ID is the same on
!!$  ! !all PEs that the corresponding berg may overlap. Then, the max/min x and y coordinates of
!!$  ! !each berg are passed between the PEs that the berg overlaps. Next, a rectangular array of
!!$  ! !bonded particles in initialized that spans these max/min coordinates, and which is defined
!!$  ! !identically and redundantly for each PE that the berg overlaps. On each PE, excess particles which
!!$  ! !do not overlap any berg cell associated with the current berg are then trimmed off. Initialization
!!$  ! !of the fully-bonded berg is completed in the iceberg module, where the bonded particles of the
!!$  ! !berg are connected across PE boundaries to produce the full iceberg.
!!$
!!$  ! !TODO: Note that the code here assumes that any two neighboring cells, each with mask>0, must
!!$  ! !be part of the same berg. If we want to calve adjacent bergs, we can calve one on the first timestep
!!$  ! !and the other on the next timestep after the first is fully initialized. Or if using damage, calve both
!!$  ! !as one berg, interp the damage to the new berg, and break bonds where there is a rift.
!!$
!!$  ! !Alternatively, if we want to allow multiple, adjacent bergs, to calve on a single
!!$  ! !time step, then we could assign different (positive, non-zero) mask values to differentiate the
!!$  ! !bergs. After the single berg is calved, we can break bonds where the mask is different (we would
!!$  ! !not break over halo cells unless the mask indicated, though the mask values may differ on other PEs).
!!$  ! !In other words, we would initialize the multiple and adjacent bergs as a single berg, which is
!!$  ! !which is subsequently broken up into the multiple adjecent bergs. This approach guarantees that
!!$  ! !the adjacent bergs are initialized without inter-berg particle overlap/separation.
!!$
!!$  ! !Alternatively, we could just have a "background" set of bonded particles that are
!!$  ! !kept constant over time, and initialize bergs by copying over a subset of these bonded
!!$  ! !particles as needed to represent the new bergs. However, this approach is not
!!$  ! !as versatile for controlling particle size. Also, it is not simple to guarantee that
!!$  ! !the particles extents of the (meridionally-aligned) first and last column would align (zonally).
!!$
!!$  ! !comment out halo-filling of tabular_calve_mask for now. If there are multiple, adjacent calving events, we
!!$  ! !may define them with different tabular_calve_mask>0 on each PE, so that the value of tabular_calve_mask may
!!$  ! !differ between PEs for the same berg. However, we make sure each PE calculates its own tabular_calve_mask
!!$  ! !in its own halo cells, and then make a copy of the resulting tabular_calve_mask. Then pass_var the original
!!$  ! !tabular_calve_mask and compare to the copy to backcalculate a consistent value for
!!$  ! !tabular_calve_mask for each berg
!!$  ! !that can be shared between each PE later.
!!$  ! !call pass_var(TC%tabular_calve_mask, G%domain)
!!$
!!$  ! !TODO: Icebergs may overlap with each other or the ice shelf, so make sure mass scaling is done appropriately
!!$  ! !on new iceberg particles to avoid issues with pressure on the ocean...Also should strongly force bergs away
!!$  ! !from ice shelves. Perhaps need to ensure that the total pressure on the ocean does not exceed that
!!$  ! !of a combination
!!$  ! !of the iceberg(s) and ice shelf mass should it fully-cover the cell (with adjustments for mass-weighting and the
!!$  ! !percentage of the total unadjusted mass that the pressure of each component exerts on the cell...). Mass will not
!!$  ! !be conserved at that instant, but ultimately is over time.
!!$
!!$  ! !1) Generate a unique label (TC%c_id) for each berg on the computational domain of a PE
!!$  ! call initialize_tabular_calving_labels_1PE(G, TC%tabular_calve_mask, TC%c_id)
!!$
!!$  ! !2) fill halos with the unique berg labels. Find connected berg cells, and update them
!!$  ! !   with the lowest berg label (TC%c_id) of all of the connected berg cells. Repeat until no changes.
!!$  ! call update_tabular_calving_labels_over_pes(G, mask, TC%c_id)
!!$
!!$  ! !3) Make a list of each of the i bergs on the local PE domain (TC%berg_list(i,1)), their
!!$  ! !   min (TC%berg_list(i,2)) max longitude (TC%berg_list(i,3)), and their
!!$  ! !   min (TC%berg_list(i,4)) max latitude (TC%berg_list(i,5))
!!$  ! call tabular_berg_info(G, TC)
!!$
!!$  ! !4) Initialize iKID icebergs over these bounds, and remove excess particles. Call this from SIS2?
!!$  ! call initialize_bonded_bergs_from_shelf(bergs, TC, h_shelf, frac_shelf_h)
!!$
!!$end subroutine process_tabular_calving
!!$
!!$!> Updates the ice shelf tabular calving mask using data from a file.
!!$subroutine update_tabular_calving_mask(G, CS, TC, Time)
!!$  type(ocean_grid_type), intent(inout) :: G   !< The ocean's grid structure.
!!$  type(ice_shelf_CS),    pointer       :: CS    !< A pointer to the control structure returned
!!$                                                 !! by a previous call to initialize_ice_shelf.
!!$  type(ice_shelf_state), intent(inout) :: TC  !< The ice shelf tabular calving state type that is being updated
!!$  type(time_type),       intent(in)    :: Time !< The current model time
!!$  ! local variables
!!$  integer :: i, j, is, ie, js, je, m, n
!!$  real, allocatable, dimension(:,:)   :: tmp2d ! Temporary array for storing ice shelf input data
!!$
!!$  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
!!$
!!$  if (TC%tabular_calving_from_file) then
!!$
!!$    if (CS%rotate_index) then
!!$      allocate(tmp2d(CS%Grid_in%isd:CS%Grid_in%ied,CS%Grid_in%jsd:CS%Grid_in%jed), source=0.0)
!!$    else
!!$      allocate(tmp2d(is:ie,js:je), source=0.0)
!!$    endif
!!$
!!$    call time_interp_external(TC%tabular_calving_handle, Time, tmp2d)
!!$    call rotate_array(tmp2d, CS%turns, TC%tabular_calve_mask)
!!$    deallocate(tmp2d)
!!$
!!$    !for now, only calve where there is ice present.
!!$    do j=js,je ; do i=is,ie
!!$      if (CS%area_shelf_h(i,j)<=0) TC%tabular_calve_mask(i,j)=0
!!$    enddo; enddo
!!$
!!$  else
!!$    !TODO: implement a tabular calving law here.
!!$    call MOM_error(FATAL, "Tabular calving is currently only possible by reading in the tabular calving mask "//&
!!$                   "from a file (see subroutine initialize_tabular_calving), i.e. there is not yet a tabular "//&
!!$                   "calving law.")
!!$  endif
!!$
!!$  call pass_var(TC%tabular_calve_mask, G%domain, complete=.true.)
!!$
!!$  !Adjust mask to make sure it extends 2 cells past the calving front. The mask will not change until the
!!$  !transition period for the gradual transition between ice shelf and iceberg is complete, so this
!!$  !extension should account for any ice that advects into these cells over this time period.
!!$  !TODO: If for some strange reason, 2 cells is not enough, you may need to just calve excess ice as
!!$  !non-interactive bergs.
!!$
!!$  !start by set calving mask to zero wherever there is no ice. This probably will not be necessary when
!!$  !a calving law is implemented, but may be helpful when reading in calving masks from file, which
!!$  !could extend into the ocean as a safeguard to make sure that the ice front is calved regardless of
!!$  !any ice front advection.
!!$  do j=js,je ; do i=is,ie
!!$    if (CS%area_shelf_h(i,j)<=0) TC%tabular_calve_mask(i,j)=0
!!$  enddo; enddo
!!$
!!$  !then set calving mask to one wherever non-ice cells are within 2 cells of a calving cell
!!$  do j=js,je; do i=is,ie
!!$    if (TC%tabular_calve_mask(i,j)>0 .and. CS%area_shelf_h(i,j)>0) then
!!$      do m=j-2,j+2; do n=i-2,i+2
!!$        if (CS%area_shelf_h(m,n)<=0) TC%tabular_calve_mask(i,j)=1.0
!!$      enddo; enddo
!!$    endif
!!$  enddo; enddo
!!$
!!$  call pass_var(TC%tabular_calve_mask, G%domain, complete=.true.)
!!$
!!$end subroutine update_tabular_calving_mask

! !> Initializes labels for the grid cells that comprise a tabular iceberg that is about to calve from an ice
! !! shelf (i.e. a group of all neighboring grid cells where calving mask > 0). Considers the current PE only.
! subroutine initialize_tabular_calving_labels_1PE(G, mask, c_id)
!   type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.
!   real, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(in)    :: mask !< A mask that is greater than zero where a
!                                                !! tabular berg should calve from a shelf
!   integer, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(in)    :: c_id !< unique ID assigned to all cells
!                                                !! of a calving tabular berg
!   integer :: i, j

!   do j=G%jsc,G%jec; do i=G%isc,G%iec
!     if (mask(i,j) .gt. 0 .and. (c_id(i,j) .eq. 0)) then
!       c_id(i,j) = generate_cell_id(G,i,j)
!       call label_tabular_bergs(G, i, j, mask, c_id)
!     endif
!   enddo; enddo
! end subroutine initialize_tabular_calving_labels_1PE

! !> Assigns the same label to all grid cells that comprise a tabular iceberg that is about to calve
! !! from an ice shelf (i.e. a group of all neighboring grid cells where calving mask > 0).
! recursive subroutine label_tabular_bergs(G, ic, jc, mask, c_id)
!   type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.
!   real, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(in)    :: mask !< A mask that is greater than zero where a
!                                                !! tabular berg should calve from a shelf
!   integer, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(in)    :: c_id !< unique ID assigned to all cells
!                                                !! of a calving tabular berg
!   integer,               intent(in)    :: ic  !< The i-index of the input cell
!   integer,               intent(in)    :: jc  !< The j-index of the input cell

!   integer :: i, j

!   do j=max(jc-1,G%jsd),min(jc+1,G%jed); do i=max(ic-1,G%isd),min(ic+1,G%ied)

!     if ((mask(i,j) > 0) .and. (c_id(i,j) .ne. c_id(ic,jc)) then

!       c_id(i,j) = min(c_id(ic,jc),c_id(i,j))
!       call label_tabular_bergs(G, i, j, mask, c_id)
!     endif
!   enddo
! end subroutine label_tabular_bergs

! !> Adjusts labels of tabular icebergs on the grid so that they are consistent between all PEs that
! !! the tabular bergs may overlap
! subroutine update_tabular_calving_labels_over_pes(G, mask, c_id)
!   type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.
!   real, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(in)    :: mask !< A mask that is greater than zero where a
!                                                !! tabular berg should calve from a shelf
!   integer, dimension(SZDI_(G),SZDJ_(G)), &
!                          intent(in)    :: c_id !< unique ID assigned to all cells
!                                                !! of a calving tabular berg
!   integer :: i, j, k, i2, j2
!   integer :: change

!   change=1
!   do while (change==1)
!     change=0
!     call pass_var(c_id, G%domain)

!     do k=1,2
!       if (k==1) then
!         i=G%isc; i2=i-1
!       else
!         i=G%iec; i2=i+1
!       endif

!       do j=G%jsc,G%jec
!         if (c_id(i,j) .ne. 0) then
!           if ((c_id(i2,j) .ne. 0) .and. (c_id(i2,j) .ne. c_id(i,j))) then
!             c_id(i,j) = min(c_id(i,j), c_id(i2,j))
!             change=1
!             call label_tabular_bergs(G, i, j, mask, c_id)
!           endif
!         endif
!       enddo

!       if (k==1) then
!         j=G%jsc; j2=j-1
!       else
!         j=G%iec; j2=j+1
!       endif

!       do i=G%isc,G%iec
!         if (c_id(i,j) .ne. 0) then
!           if ((c_id(i,j2) .ne. 0) .and. (c_id(i,j2) .ne. c_id(i,j))) then
!             c_id(i,j) = min(c_id(i,j), c_id(i,j2))
!             change=1
!             call label_tabular_bergs(G, i, j, mask, c_id)
!           endif
!         endif
!       enddo
!     enddo
!     call max_across_PEs(change)
!   enddo

! end subroutine update_tabular_calving_labels_over_pes

! !> Finds the unique bergs of a computational domain field, and their global extent
! !! (min and max latitude and longitude)
! subroutine tabular_berg_info(G, TC)
!   type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
!   type(tabular_calving_state), pointer, intent(inout) :: TC !< A pointer to the tabular calving structure
!   integer, intent(out) :: bcount !< number of unique bergs
!   ! local variables
!   real, dimension(:), allocatable :: tmp
!   integer, pointer :: c_id(:,:)
!   integer :: n, i
!   real :: min_val, max_val
!   integer :: pe_N,pe_S,pe_E,pe_W

!   c_id=>TC%cid
!   n=(G%iec-G%isc) * (G%jec-G%jsc)
!   allocate(tmp(n))
!   tmp(:) = 0.
!   min_val = minval(c_id(G%isc:G%iec,G%jsc:G%jec))-1
!   max_val = maxval(c_id(G%isc:G%iec,G%jsc:G%jec))
!   bcount = 0
!   do while (min_val<max_val)
!     bcount = bcount+1
!     min_val = minval(c_id(G%isc:G%iec,G%jsc:G%jec), mask=c_id(G%isc:G%iec,G%jsc:G%jec)>min_val)
!     tmp(bcount) = real(min_val)
!   enddo

!   if (allocated(TC%berg_list)) deallocate(TC%berg_list)

!   allocate(TC%berg_list(bcount,5))

!   do i = 1,bcount

!     !the unique berg
!     TC%berg_list(i,1) = float(tmp(i))

!     !minlon
!     TC%berg_list(i,2) = minval(G%geoloncu(G%isc-1:G%iec-1,G%jsc:G%jec), &
!       mask=c_id(G%isc:G%iec,G%jsc:G%jec)==tmp(i))
!     !maxlon
!     TC%berg_list(i,3) = maxval(G%geoloncu(G%isc:G%iec,G%jsc:G%jec), &
!       mask=c_id(G%isc:G%iec,G%jsc:G%jec)==tmp(i))
!     !minlat
!     TC%berg_list(i,4) = minval(G%geolatcu(G%isc:G%iec,G%jsc-1:G%jec-1), &
!       mask=c_id(G%isc:G%iec,G%jsc:G%jec)==tmp(i))
!     !maxlat
!     TC%berg_list(i,5) = maxval(G%geolatcu(G%isc:G%iec,G%jsc:G%jec), &
!       mask=c_id(G%isc:G%iec,G%jsc:G%jec)==tmp(i))
!   enddo

!   ! call mpp_get_neighbor_pe(grd%domain, NORTH, grd%pe_N)
!   ! call mpp_get_neighbor_pe(grd%domain, SOUTH, grd%pe_S)
!   ! call mpp_get_neighbor_pe(grd%domain, EAST, grd%pe_E)
!   ! call mpp_get_neighbor_pe(grd%domain, WEST, grd%pe_W)

!   TC%berg_pe_count=bcount

!   !update the local coordinate bounds for each berg in TC%berg_list with the
!   !global coordinate bounds for each berg across all PEs
!   call update_berg_lists_on_all_pes(G, TC)

!   deallocate(tmp)
! end subroutine tabular_berg_info

! !> Updates TC%berg_list
! recursive subroutine update_berg_lists_on_all_pes(G, TC)
!   type(ocean_grid_type), intent(in) :: G    !< The ocean's grid structure.
!   type(tabular_calving_state), pointer, intent(inout) :: TC !< A pointer to the tabular calving structure
!   integer :: i1(4), i2(4), j1(4), j2(4), nbergs(4)
!   integer, allocatable :: tracker(:,:)
!   real, allocatable :: buffer(:)
!   integer :: n,k
!   integer :: changes, localchanges, nbergs_rcvd
!   real, pointer :: c_id(:,:), pebl(:,:)


!   !Halo range for each cardinal direction

!   i1(1)=G%iec+1; i2(1)=G%ied;   j1(1)=G%jsd;   j2(1)=G%jed   !E
!   i1(2)=G%isd;   i2(2)=G%isc-1; j1(2)=G%jsd;   j2(2)=G%jed   !W
!   i1(3)=G%isd;   i2(3)=G%ied;   j1(3)=G%jec+1; j2(3)=G%jed   !N
!   i1(4)=G%isd;   i2(4)=G%ied;   j1(4)=G%jsd;   j2(4)=G%jsc-1 !S

!   c_id=>TC%cid
!   pebl=>TC%berg_list

!   !each row of tracker is a berg on the PE, and each column corresponds to a cardinal direction.
!   !where tracker == 1, the berg in that row is to be sent to the PE in the corresponding direction.
!   allocate(tracker(TC%berg_pe_count,4); tracker=0

!   do n = 1,TC%berg_pe_count
!     do k=1,4
!       if (any(c_id(i1(k):i2(k),j1(k):j2(k))==pebl(n,1))) tracker(n,k)=1
!       if (any(c_id(i1(k):i2(k),j1(k):j2(k))==pebl(n,1))) tracker(n,k)=1
!     enddo
!   enddo

!   !number of bergs to send to each direction
!   nbergs(1:4)=sum(tracker,2)

!   changes = 1
!   localchanges = 1
!   do while (changes.ne.0)

!     changes=0

!     !send bergs east/west
!     if (G%pe_E.ne.NULL_PE) then
!       if (localchanges>0) then
!         call mpp_send(nbergs(1)*5, plen=1, to_pe=G%pe_E, tag=COMM_TAG_1)
!         if (nbergs(1).gt.0) then
!           allocate(buffer(nbergs(1)*5)
!           call pack_tabular_buffer(TC%berg_pe_count, nbergs(1)*5, 1, pebl, tracker, buffer)
!           call mpp_send(buffer, nbergs(1)*5, G%pe_E, tag=COMM_TAG_2)
!           deallocate(buffer)
!         endif
!       else
!         call mpp_send(0, plen=1, to_pe=G%pe_E, tag=COMM_TAG_1)
!       endif
!     endif
!     if (G%pe_W.ne.NULL_PE) then
!       if (localchanges>0) then
!         call mpp_send(nbergs(2)*5, plen=1, to_pe=G%pe_W, tag=COMM_TAG_1)
!         if (nbergs(2).gt.0) then
!           allocate(buffer(nbergs(2)*5)
!           call pack_tabular_buffer(TC%berg_pe_count, nbergs(2)*5, 2, pebl, tracker, buffer)
!           call mpp_send(buffer, nbergs(2)*5, G%pe_W, tag=COMM_TAG_2)
!           deallocate(buffer)
!         endif
!       else
!         call mpp_send(0, plen=1, to_pe=G%pe_W, tag=COMM_TAG_1)
!       endif
!     endif

!     !receive bergs from west/east
!     if (grd%pe_W.ne.NULL_PE) then
!       data_rcvd=0
!       call mpp_recv(data_rcvd, glen=1, from_pe=grd%pe_W, tag=COMM_TAG_1)
!       if (data_rcvd.gt.0) then
!         allocate(buffer(data_rcvd))
!         call mpp_recv(buffer, data_rcvd, grd%pe_W, tag=COMM_TAG_2)
!         call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd, pebl, buffer. changes)
!         deallocate(buffer)
!       endif
!     endif
!     if (grd%pe_E.ne.NULL_PE) then
!       data_rcvd=0
!       call mpp_recv(data_rcvd, glen=1, from_pe=grd%pe_E, tag=COMM_TAG_1)
!       if (data_rcvd.gt.0) then
!         allocate(buffer(data_rcvd))
!         call mpp_recv(buffer, data_rcvd, grd%pe_E, tag=COMM_TAG_2)
!         call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd, pebl, buffer. changes)
!         deallocate(buffer)
!       endif
!     endif

!     localchanges=max(localchanges,changes)

!     !send bergs north/south
!     if (G%pe_N.ne.NULL_PE) then
!       if (localchanges>0) then
!         call mpp_send(nbergs(3)*5, plen=1, to_pe=G%pe_N, tag=COMM_TAG_1)
!         if (nbergs(3).gt.0) then
!           allocate(buffer(nbergs(3)*5)
!           call pack_tabular_buffer(TC%berg_pe_count, nbergs(3)*5, 3, pebl, tracker, buffer)
!           call mpp_send(buffer, nbergs(3)*5, G%pe_N, tag=COMM_TAG_2)
!           deallocate(buffer)
!         endif
!       endif
!     else
!       call mpp_send(0, plen=1, to_pe=G%pe_N, tag=COMM_TAG_1)
!     endif
!     if (G%pe_S.ne.NULL_PE) then
!       if (localchanges>0) then
!         call mpp_send(nbergs(4)*5, plen=1, to_pe=G%pe_S, tag=COMM_TAG_1)
!         if (nbergs(4).gt.0) then
!           allocate(buffer(nbergs(4)*5)
!           call pack_tabular_buffer(TC%berg_pe_count, nbergs(4)*5, 4, pebl, tracker, buffer)
!           call mpp_send(buffer, nbergs(4)*5, G%pe_S, tag=COMM_TAG_2)
!           deallocate(buffer)
!         endif
!       endif
!     else
!       call mpp_send(0, plen=1, to_pe=G%pe_S, tag=COMM_TAG_1)
!     endif

!     !receive bergs north/south
!     if (grd%pe_S.ne.NULL_PE) then
!       data_rcvd=0
!       call mpp_recv(data_rcvd, glen=1, from_pe=grd%pe_S, tag=COMM_TAG_1)
!       if (data_rcvd.gt.0) then
!         allocate(buffer(data_rcvd))
!         call mpp_recv(buffer, data_rcvd, grd%pe_S, tag=COMM_TAG_2)
!         call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd, pebl, buffer. changes)
!         deallocate(buffer)
!       endif
!     endif
!     if (grd%pe_N.ne.NULL_PE) then
!       data_rcvd=0
!       call mpp_recv(data_rcvd, glen=1, from_pe=grd%pe_N, tag=COMM_TAG_1)
!       if (data_rcvd.gt.0) then
!         allocate(buffer(data_rcvd))
!         call mpp_recv(buffer, data_rcvd, grd%pe_N, tag=COMM_TAG_2)
!         call unpack_tabular_buffer_and_update_bounds(TC%berg_pe_count, data_rcvd, pebl, buffer. changes)
!         deallocate(buffer)
!       endif
!     endif

!     localchanges=changes
!     call max_across_PEs(changes)
!   enddo

! end subroutine update_berg_lists_on_all_pes

! !> pack the buffer with the info for the bergs in the c_id array being sent to another PE
! subroutine pack_tabular_buffer(lbergs, bberg_data, dir, pebl, tracker, buffer)
!   integer :: lbergs !< number of local bergs on the current PE
!   integer :: bberg_data !< count of berg data being sent in the buffer to another PE
!   integer :: dir !< direction the bergs are being sent
!   real, pointer :: pebl(:,:) !< array of bergs on the current PE and their lat/lon bounds
!   integer :: tracker(totbergs,4) !< tracks the direction to send bergs from the current PE
!   real :: buffer(bberg_data) !< the buffer being packed with bergs to send to another PE
!   integer :: k, i

!   i=1
!   do k=1,lbergs
!     if (tracker(k,dir)==1) then
!       buffer(i:i+4)=pebl(k,1:5)
!       i=i+5
!     endif
!   enddo

! end subroutine pack_tabular_buffer

! !> Unpack the buffer with the info for the bergs in the c_id array being sent from another PE
! !! if unpacked berg has more extreme bounds than the same berg in the current PE, update the
! !! bounds of the current PE berg to match
! subroutine unpack_tabular_buffer_and_update_bounds(lbergs, bberg_data, pebl, buffer, changes)
!   integer :: lbergs !< number of local bergs on the current PE
!   integer :: bberg_data !< count of berg data being received in the buffer from another PE
!   real, pointer :: pebl !< array of bergs on the current PE and their lat/lon bounds
!   real :: buffer(bberg_data) !< the buffer of bergs being received on the current PE
!   real :: buffer2(bberg_data/5,5)
!   integer :: m, n, changes

!   !reshape the buffer
!   do n=1,bberg_data/5
!     buffer2(n,1:5) = buffer((n-1)*5+1:(n-1)*5+5)
!   enddo

!   do m=1,lbergs
!     do n=1,bbergs
!       if (pebl(m,1)==buffer2(n,1)) then
!         !same berg ID detected in current PE berg list and berg buffer from other PE
!         !extend current berg bounds if needed, to reflect more extensive bounds from other PE
!         if (pebl(m,2)>buffer2(n,2)) then; pebl(m,2)=buffer2(n,2); changes=changes+1; endif
!         if (pebl(m,3)<buffer2(n,3)) then; pebl(m,3)=buffer2(n,3); changes=changes+1; endif
!         if (pebl(m,4)>buffer2(n,4)) then; pebl(m,4)=buffer2(n,4); changes=changes+1; endif
!         if (pebl(m,5)<buffer2(n,5)) then; pebl(m,5)=buffer2(n,5); changes=changes+1; endif
!       endif
!     enddo
!   enddo

! end subroutine unpack_tabular_buffer_and_update_bounds

! !> Calculate cell_id, a unique integer for each grid cell.
! !! This is the same as ij_component_of_id in icebergs_framework.
! xinteger function generate_cell_id(G, i, j)
!   type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.
!   integer,                intent(in) :: i   !< i-index of grid cell
!   integer,                intent(in) :: j   !< j-index of grid cell
!   ! Local variables
!   integer :: iNg ! Zonal size of the global grid

!   ! Using the current grid shape maximizes the numbers of IDs that can be represented
!   ! allowing up to 30-minute uniform global resolution, or potentially finer if non-uniform.
!   iNg = G%ieg - G%isg + 1

!   ! ij_component_of_id is unique number for each grid cell
!   !(32-bit integers allow for ~1/100th degree global resolution)
!   generate_cell_id = i + ( iNg * ( j - 1 ) )
! end function generate_cell_id

end module MOM_ice_shelf_tabular_calving
