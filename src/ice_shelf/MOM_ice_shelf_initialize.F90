! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initialize ice shelf variables
module MOM_ice_shelf_initialize

use MOM_grid, only : ocean_grid_type
use MOM_array_transform,      only : rotate_array
use MOM_hor_index,  only : hor_index_type
use MOM_file_parser, only : get_param, read_param, log_param, param_file_type
use MOM_io, only: MOM_read_data, file_exists, field_exists, slasher, CORNER
use MOM_error_handler, only : MOM_error, MOM_mesg, FATAL, WARNING, is_root_pe
use MOM_unit_scaling, only : unit_scale_type
use user_shelf_init, only: USER_init_ice_thickness
use MOM_domains, only : pass_var, CORNER

implicit none ; private

#include <MOM_memory.h>

public initialize_ice_thickness
public initialize_ice_shelf_boundary_channel
public initialize_ice_flow_from_file
public initialize_bed_node_from_file
public corner_cell_weights, nodal_cell_mean
public initialize_ice_shelf_boundary_from_file
public initialize_ice_C_basal_friction
public initialize_ice_AGlen
public initialize_ice_SMB
! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Initialize ice shelf thickness
subroutine initialize_ice_thickness(h_shelf, area_shelf_h, hmask, melt_mask, G, G_in, US, PF, &
                                    rotate_index, turns, h_nodal)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  type(ocean_grid_type), intent(in)    :: G_in    !< The ocean's unrotated grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: melt_mask !< A mask indicating where to allow ice-shelf melting [nondim]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters
  logical, intent(in), optional        :: rotate_index !< If true, this is a rotation test
  integer, intent(in), optional        :: turns !< Number of turns for rotation test
  real, dimension(:,:,:,:), optional, pointer :: h_nodal !< DG(1) corner thickness, filled only
                                             !! under INIT_ICE_THICKNESS_NODAL [Z ~> m].

  character(len=40)  :: mdl = "initialize_ice_thickness" ! This subroutine's name.
  character(len=200) :: config
  logical :: rotate = .false.
  logical :: nodal_thickness ! If true, the thickness is read at B-grid corner nodes.
  real, allocatable, dimension(:,:) :: tmp1_2d ! Temporary array for storing ice shelf input data [Z~>m]
  real, allocatable, dimension(:,:) :: tmp2_2d ! Temporary array for storing ice shelf input data [L2~>m2]
  real, allocatable, dimension(:,:) :: tmp3_2d ! Temporary array for storing ice shelf input data [nondim]
  real, allocatable, dimension(:,:) :: tmp4_2d ! Temporary array for storing ice shelf input data [nondim]

  call get_param(PF, mdl, "ICE_PROFILE_CONFIG", config, &
                 "This specifies how the initial ice profile is specified. "//&
                 "Valid values are: CHANNEL, FILE, and USER.", &
                 fail_if_missing=.true.)

  if (PRESENT(rotate_index)) rotate=rotate_index

  if (rotate) then
    ! The nodal thickness would need a corner-position rotate_array, which does not exist.
    call get_param(PF, mdl, "INIT_ICE_THICKNESS_NODAL", nodal_thickness, &
                   default=.false., do_not_log=.true.)
    if (nodal_thickness) call MOM_error(FATAL, "initialize_ice_thickness: "//&
                 "INIT_ICE_THICKNESS_NODAL is not supported in a rotation test.")
    allocate(tmp1_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=0.0)
    allocate(tmp2_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=0.0)
    allocate(tmp3_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=0.0)
    allocate(tmp4_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=1.0)
    select case ( trim(config) )
      case ("CHANNEL") ; call initialize_ice_thickness_channel (tmp1_2d, tmp2_2d, tmp3_2d, G_in, US, PF)
      case ("FILE") ; call initialize_ice_thickness_from_file (tmp1_2d, tmp2_2d, tmp3_2d, tmp4_2d, &
                                                               G_in, US, PF)
      case ("USER") ; call USER_init_ice_thickness (tmp1_2d, tmp2_2d, tmp3_2d, G_in, US, PF)
      case default  ; call MOM_error(FATAL,"MOM_initialize: Unrecognized ice profile setup "//trim(config))
    end select
    call rotate_array(tmp1_2d,turns, h_shelf)
    call rotate_array(tmp2_2d,turns, area_shelf_h)
    call rotate_array(tmp3_2d,turns, hmask)
    call rotate_array(tmp4_2d,turns, melt_mask)
    deallocate(tmp1_2d,tmp2_2d,tmp3_2d,tmp4_2d)
  else
    select case ( trim(config) )
      case ("CHANNEL") ; call initialize_ice_thickness_channel (h_shelf, area_shelf_h, hmask, G, US, PF)
      case ("FILE") ; call initialize_ice_thickness_from_file (h_shelf, area_shelf_h, hmask, melt_mask, &
                                                               G, US, PF, h_nodal)
      case ("USER") ; call USER_init_ice_thickness (h_shelf, area_shelf_h, hmask, G, US, PF)
      case default  ; call MOM_error(FATAL,"MOM_initialize: Unrecognized ice profile setup "//trim(config))
    end select
  endif

end subroutine initialize_ice_thickness

!> Initialize ice shelf thickness from file
subroutine initialize_ice_thickness_from_file(h_shelf, area_shelf_h, hmask, melt_mask, G, US, PF, h_nodal)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: melt_mask !< A mask indicating where to allow ice-shelf melting [nondim]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters
  real, dimension(:,:,:,:), optional, pointer :: h_nodal !< DG(1) corner thickness, filled only
                                             !! under INIT_ICE_THICKNESS_NODAL [Z ~> m].

  !  This subroutine reads ice thickness and area from a file and puts it into
  !  h_shelf [Z ~> m] and area_shelf_h [L2 ~> m2] (and dimensionless) and updates hmask
  character(len=200) :: filename,thickness_file,inputdir ! Strings for file/path
  character(len=200) :: thickness_varname, area_varname, hmask_varname, melt_mask_varname  ! Variable name in file
  character(len=40)  :: mdl = "initialize_ice_thickness_from_file" ! This subroutine's name.
  real, allocatable :: h_node(:,:) ! Nodal (B-grid corner) ice thickness [Z ~> m].
  integer :: i, j, isc, jsc, iec, jec
  logical :: hmask_set
  logical :: nodal_thickness ! If true, read the thickness at B-grid corner nodes.
  logical :: reentrant_x, reentrant_y ! True if the domain wraps in x or y
  real :: len_sidestress, udh

  call MOM_mesg("Initialize_ice_thickness_from_file: reading thickness")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "ICE_THICKNESS_FILE", thickness_file, &
                 "The file from which the bathymetry is read.", &
                 default="ice_shelf_h.nc")
  call get_param(PF, mdl, "LEN_SIDE_STRESS", len_sidestress, &
                 "position past which shelf sides are stress free.", &
                 default=0.0, units="axis_units")

  filename = trim(inputdir)//trim(thickness_file)
  call log_param(PF, mdl, "INPUTDIR/THICKNESS_FILE", filename)
  call get_param(PF, mdl, "ICE_THICKNESS_VARNAME", thickness_varname, &
                 "The name of the thickness variable in ICE_THICKNESS_FILE. With "//&
                 "INIT_ICE_THICKNESS_NODAL this names a B-grid corner field instead of "//&
                 "a cell-averaged one.", &
                 default="h_shelf")
  call get_param(PF, mdl, "INIT_ICE_THICKNESS_NODAL", nodal_thickness, &
                 "If true, read ICE_THICKNESS_VARNAME at B-grid nodes and set the cell thickness "//&
                 "to the area-weighted mean of their bilinear interpolant. The nodes also "//&
                 "initialize USE_DG_THICKNESS. Requires the shelf mask in the file.", &
                 default=.false.)
  call get_param(PF, mdl, "ICE_AREA_VARNAME", area_varname, &
                 "The name of the area variable in ICE_THICKNESS_FILE.", &
                 default="area_shelf_h")
  hmask_varname="h_mask"
  call get_param(PF, mdl, "MELT_MASK_VARNAME", melt_mask_varname, &
                 "The name of the melt mask variable in ICE_THICKNESS_FILE.", &
                 default="melt_mask")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))
  if (.not.nodal_thickness) &
    call MOM_read_data(filename, trim(thickness_varname), h_shelf, G%Domain, scale=US%m_to_Z)
  call MOM_read_data(filename,trim(area_varname), area_shelf_h, G%Domain, scale=US%m_to_L**2)
  if (field_exists(filename, trim(hmask_varname), MOM_domain=G%Domain)) then
    call MOM_read_data(filename, trim(hmask_varname), hmask, G%Domain)
    hmask_set = .true.
  else
    call MOM_error(WARNING, "Ice shelf thickness initialized without setting the shelf mask "//&
              "from variable "//trim(hmask_varname)//", which does not exist in "//trim(filename))
    hmask_set = .false.
  endif
  if (field_exists(filename, trim(melt_mask_varname), MOM_domain=G%Domain)) then
    call MOM_read_data(filename, trim(melt_mask_varname), melt_mask, G%Domain)
  else
    melt_mask(:,:)=1.0
  endif

  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  if (.not.hmask_set) then
    ! Set hmask based on the values in h_shelf.
    if (nodal_thickness) call MOM_error(FATAL, "initialize_ice_thickness_from_file: "//&
              "INIT_ICE_THICKNESS_NODAL derives hmask from no cell-averaged thickness, so "//&
              "it requires the variable "//trim(hmask_varname)//" in "//trim(filename)//".")
    do j=jsc,jec ; do i=isc,iec
      hmask(i,j) = 0.0
      if (h_shelf(i,j) > 0.0) hmask(i,j) = 1.0
    enddo ; enddo
  endif

  ! Nodal thickness and its cell means, once hmask is known.
  if (nodal_thickness) then
    if (.not. present(h_nodal)) call MOM_error(FATAL, "initialize_ice_thickness_from_file: "//&
              "INIT_ICE_THICKNESS_NODAL needs the nodal thickness array, which only the "//&
              "dynamic ice-shelf initialization path provides.")
    if (.not. associated(h_nodal)) call MOM_error(FATAL, "initialize_ice_thickness_from_file: "//&
              "INIT_ICE_THICKNESS_NODAL was requested before the nodal thickness was allocated.")
    if (len_sidestress > 0.) call MOM_error(FATAL, "initialize_ice_thickness_from_file: "//&
              "LEN_SIDE_STRESS tapers the cell-averaged thickness, which "//&
              "INIT_ICE_THICKNESS_NODAL derives from the nodal field rather than reading, "//&
              "leaving the taper and the nodal thickness inconsistent.")
    call get_param(PF, mdl, "REENTRANT_X", reentrant_x, &
                 " If true, the domain is zonally reentrant.", &
                 default=.false., do_not_log=.true.)
    call get_param(PF, mdl, "REENTRANT_Y", reentrant_y, &
                 " If true, the domain is meridionally reentrant.", &
                 default=.false., do_not_log=.true.)
    allocate(h_node(G%IsdB:G%IedB, G%JsdB:G%JedB), source=0.0)
    call MOM_read_data(filename, trim(thickness_varname), h_node, G%Domain, &
                       position=CORNER, scale=US%m_to_Z)
    call pass_var(h_node, G%domain, position=CORNER)

    h_nodal(:,:,:,:) = 0.0
    do j=jsc,jec ; do i=isc,iec
      if (hmask(i,j)==0) then
        h_shelf(i,j) = 0.0
      else
        h_nodal(i,j,1,1) = h_node(I-1, J-1)
        h_nodal(i,j,2,1) = h_node(I,   J-1)
        h_nodal(i,j,1,2) = h_node(I-1, J  )
        h_nodal(i,j,2,2) = h_node(I,   J  )
        h_shelf(i,j) = corner_cell_mean(G, i, j, reentrant_x, reentrant_y, &
                                        h_nodal(i,j,1,1), h_nodal(i,j,2,1), &
                                        h_nodal(i,j,1,2), h_nodal(i,j,2,2))
      endif
    enddo ; enddo
    deallocate(h_node)
  endif

    do j=jsc,jec
      do i=isc,iec

      ! taper ice shelf in area where there is no sidestress -
      ! but do not interfere with hmask

        if ((len_sidestress > 0.) .and. (G%geoLonCv(i,j) > len_sidestress)) then
          udh = exp(-(G%geoLonCv(i,j)-len_sidestress)/5.0) * h_shelf(i,j)
          if (udh <= 25.0) then
            h_shelf(i,j) = 0.0
            area_shelf_h(i,j) = 0.0
          else
            h_shelf(i,j) = udh
          endif
        endif

        ! Fix for any round-off difference in ice-shelf area between the file and model grid
        if (hmask_set) then
          if (hmask(i,j)==1 .or. hmask(i,j)==3) area_shelf_h(i,j)=G%areaT(i,j)
        endif

      ! update thickness mask

        if (area_shelf_h(i,j) >= G%areaT(i,j)) then
          if (.not. hmask_set) hmask(i,j) = 1.
          area_shelf_h(i,j)=G%areaT(i,j)
        elseif (area_shelf_h(i,j) == 0.0) then
          hmask(i,j) = 0.
        elseif ((area_shelf_h(i,j) > 0) .and. (area_shelf_h(i,j) <= G%areaT(i,j))) then
          hmask(i,j) = 2.
        else
          call MOM_error(FATAL,mdl// " AREA IN CELL OUT OF RANGE")
        endif
      enddo
    enddo
end subroutine initialize_ice_thickness_from_file

!> Initialize ice shelf thickness for a channel configuration
subroutine initialize_ice_thickness_channel(h_shelf, area_shelf_h, hmask, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< The ice shelf thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: area_shelf_h !< The area per cell covered by the ice shelf [L2 ~> m2].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  character(len=40)  :: mdl = "initialize_ice_shelf_thickness_channel" ! This subroutine's name.
  real :: max_draft, min_draft, flat_shelf_width, c1, slope_pos
  real :: edge_pos, shelf_slope_scale
  integer :: i, j, jsc, jec, jsd, jed, jedg, nyh, isc, iec, isd, ied
  integer :: j_off

  jsc = G%jsc ; jec = G%jec ; isc = G%isc ; iec = G%iec
  jsd = G%jsd ; jed = G%jed ; isd = G%isd ; ied = G%ied
  nyh = G%domain%njhalo ; jedg = G%domain%njglobal+nyh
  j_off = G%jdg_offset

  call MOM_mesg(mdl//": setting thickness")

  call get_param(PF, mdl, "SHELF_MAX_DRAFT", max_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(PF, mdl, "SHELF_MIN_DRAFT", min_draft, &
                 units="m", default=1.0, scale=US%m_to_Z)
  call get_param(PF, mdl, "FLAT_SHELF_WIDTH", flat_shelf_width, &
                 units="axis_units", default=0.0)
  call get_param(PF, mdl, "SHELF_SLOPE_SCALE", shelf_slope_scale, &
                 units="axis_units", default=0.0)
  call get_param(PF, mdl, "SHELF_EDGE_POS_0", edge_pos, &
                 units="axis_units", default=0.0)
!  call get_param(param_file, mdl, "RHO_0", Rho_ocean, &
!                 "The mean ocean density used with BOUSSINESQ true to "//&
!                 "calculate accelerations and the mass for conservation "//&
!                 "properties, or with BOUSSINESQ false to convert some "//&
!                 "parameters from vertical units of m to kg m-2.", &
!                 units="kg m-3", default=1035.0, scale=US%Z_to_m)

  slope_pos = edge_pos - flat_shelf_width
  c1 = 0.0 ; if (shelf_slope_scale > 0.0) c1 = 1.0 / shelf_slope_scale


  do j=G%jsd,G%jed

  if (((j+j_off) <= jedg) .AND. ((j+j_off) >= nyh+1)) then

    do i=G%isc,G%iec

      if ((j >= jsc) .and. (j <= jec)) then

        if (G%geoLonCu(i-1,j) >= edge_pos) then
        ! Everything past the edge is open ocean.
          area_shelf_h(i,j) = 0.0
          hmask (i,j) = 0.0
          h_shelf (i,j) = 0.0
        else
          if (G%geoLonCu(i,j) > edge_pos) then
            area_shelf_h(i,j) = G%areaT(i,j) * (edge_pos - G%geoLonCu(i-1,j)) / &
                                (G%geoLonCu(i,j) - G%geoLonCu(i-1,j))
            hmask (i,j) = 2.0
          else
            area_shelf_h(i,j) = G%areaT(i,j)
            hmask (i,j) = 1.0
          endif

          if (G%geoLonT(i,j) > slope_pos) then
            h_shelf(i,j) = min_draft
          else
            h_shelf(i,j) = (min_draft + &
               (max_draft - min_draft) * &
               min(1.0, (c1*(slope_pos - G%geoLonT(i,j)))**2) )
          endif

        endif
      endif

      if ((i+G%idg_offset) == G%domain%nihalo+1) then
        hmask(i-1,j) = 3.0
      endif

    enddo
  endif ; enddo

end subroutine initialize_ice_thickness_channel

!> Initialize ice shelf boundary conditions for a channel configuration
subroutine initialize_ice_shelf_boundary_channel(u_face_mask_bdry, v_face_mask_bdry, &
                u_flux_bdry_val, v_flux_bdry_val, u_bdry_val, v_bdry_val, u_shelf, v_shelf, h_bdry_val, &
                hmask,  h_shelf, G, US, PF )

  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_face_mask_bdry !< A boundary-type mask at C-grid u faces

  real, dimension(SZIB_(G),SZJ_(G)), &
                         intent(inout) :: u_flux_bdry_val  !< The boundary thickness flux through
                                                     !! C-grid u faces [L Z T-1 ~> m2 s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_face_mask_bdry !< A boundary-type mask at C-grid v faces

  real, dimension(SZI_(G),SZJB_(G)), &
                         intent(inout) :: v_flux_bdry_val  !< The boundary thickness flux through
                                                     !! C-grid v faces [L Z T-1 ~> m2 s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_bdry_val !< The zonal ice shelf velocity at open
                                                      !! boundary vertices [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_bdry_val !< The meridional ice shelf velocity at open
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_shelf !< The zonal ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_shelf !< The meridional ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_bdry_val !< The ice shelf thickness at open boundaries [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_shelf !< Ice-shelf thickness [Z ~> m]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  character(len=40)  :: mdl = "initialize_ice_shelf_boundary_channel" ! This subroutine's name.
  integer :: i, j, isd, jsd, giec, gjec, gisc, gjsc,gisd,gjsd, isc, jsc, iec, jec, ied, jed
  real    :: input_thick ! The input ice shelf thickness [Z ~> m]
  real    :: input_vel  ! The input ice velocity at the upstream boundary [L T-1 ~> m s-1]
  real    :: lenlat, len_stress, westlon, lenlon, southlat ! The input positions of the channel boundarises

  lenlat = G%len_lat
  lenlon = G%len_lon
  westlon = G%west_lon
  southlat = G%south_lat

  call get_param(PF, mdl, "INPUT_VEL_ICE_SHELF", input_vel, &
                 "inflow ice velocity at upstream boundary", &
                 units="m s-1", default=0., scale=US%m_s_to_L_T)
  call get_param(PF, mdl, "INPUT_THICK_ICE_SHELF", input_thick, &
                 "flux thickness at upstream boundary", &
                 units="m", default=1000., scale=US%m_to_Z)
  call get_param(PF, mdl, "LEN_SIDE_STRESS", len_stress, &
                 "maximum position of no-flow condition in along-flow direction", &
                 units="km", default=0.)

  call MOM_mesg(mdl//": setting boundary")

  isd = G%isd ; ied = G%ied
  jsd = G%jsd ; jed = G%jed
  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec
  gjsd = G%Domain%njglobal ; gisd = G%Domain%niglobal
  gisc = G%Domain%nihalo ; gjsc = G%Domain%njhalo
  giec = G%Domain%niglobal+gisc ; gjec = G%Domain%njglobal+gjsc

 !---------b.c.s based on geopositions -----------------
  do j=jsc,jec+1
    do i=isc-1,iec+1
 ! upstream boundary - set either dirichlet or flux condition

      if (G%geoLonBu(i,j) == westlon) then
        hmask(i+1,j) = 3.0
        !---
        !OLD: thickness_bdry_val was used for ice dynamics, and h_bdry_val was not used anywhere except here:
        !h_bdry_val(i+1,j) = h_shelf(i+1,j) ; thickness_bdry_val(i+1,j) = h_bdry_val(i+0*1,j)
        !---
        !NEW: h_bdry_val is used for ice dynamics instead of thickness_bdry_val, which was removed
        h_bdry_val(i+1,j) = h_shelf(i+0*1,j) !why 0*1
        !---
        u_face_mask_bdry(i+1,j) = 5.0
        u_bdry_val(i+1,j) = input_vel*(1-16.0*((G%geoLatBu(i-1,j)/lenlat-0.5))**4) !velocity distribution
      endif


      ! side boundaries: no flow
      if (G%geoLatBu(i,j-1) == southlat) then !bot boundary
        if (len_stress == 0. .OR. G%geoLonCv(i,j) <= len_stress) then
          v_face_mask_bdry(i,j+1) = 0.
          u_face_mask_bdry(i,j) = 3.
          u_bdry_val(i,j) = 0.
          v_bdry_val(i,j) = 0.
        else
          v_face_mask_bdry(i,j+1) = 1.
          u_face_mask_bdry(i,j) = 3.
          u_bdry_val(i,j) = 0.
          v_bdry_val(i,j) = 0.
        endif
      elseif (G%geoLatBu(i,j-1) == southlat+lenlat) then !top boundary
        if (len_stress == 0. .OR. G%geoLonCv(i,j) <= len_stress) then
          v_face_mask_bdry(i,j-1) = 0.
          u_face_mask_bdry(i,j-1) = 3.
        else
          v_face_mask_bdry(i,j-1) = 3.
          u_face_mask_bdry(i,j-1) = 3.
        endif
      endif

      ! downstream boundary - CFBC
      if (G%geoLonBu(i,j) == westlon+lenlon) then
        u_face_mask_bdry(i-1,j) = 2.0
      endif

    enddo
  enddo
end subroutine initialize_ice_shelf_boundary_channel


!> Initialize ice shelf flow from file
subroutine initialize_ice_flow_from_file(bed_elev,u_shelf, v_shelf,float_cond,&
                                         G, US, PF, skip_bed)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: bed_elev !< The bed elevation   [Z ~> m].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_shelf !< The zonal ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_shelf !< The meridional ice shelf velocity  [L T-1 ~> m s-1].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout)    :: float_cond !< An array indicating where the ice
                                                !! shelf is floating: 0 if floating, 1 if not. [nondim]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters
  logical, optional,     intent(in)    :: skip_bed !< If true, do not read the cell bed,
                                                !! because it comes from the nodal bed.

  logical :: skip_bed_local

  !  This subroutine reads ice thickness and area from a file and puts it into
  !  h_shelf [Z ~> m] and area_shelf_h [L2 ~> m2] (and dimensionless) and updates hmask
  character(len=200) :: filename,vel_file,inputdir,bed_topo_file ! Strings for file/path
  character(len=200) :: ushelf_varname, vshelf_varname, &
                        floatfr_varname, bed_varname  ! Variable name in file
  character(len=40)  :: mdl = "initialize_ice_velocity_from_file" ! This subroutine's name.

  call MOM_mesg("  MOM_ice_shelf_init_profile.F90, initialize_velocity_from_file: reading velocity")

  skip_bed_local = .false.
  if (present(skip_bed)) skip_bed_local = skip_bed

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "ICE_VELOCITY_FILE", vel_file, &
                 "The file from which the velocity is read.", &
                 default="ice_shelf_vel.nc")

  filename = trim(inputdir)//trim(vel_file)
  call log_param(PF, mdl, "INPUTDIR/THICKNESS_FILE", filename)
  call get_param(PF, mdl, "ICE_U_VEL_VARNAME", ushelf_varname, &
                 "The name of the u velocity variable in ICE_VELOCITY_FILE.", &
                 default="u_shelf")
  call get_param(PF, mdl, "ICE_V_VEL_VARNAME", vshelf_varname, &
                 "The name of the v velocity variable in ICE_VELOCITY_FILE.", &
                 default="v_shelf")
  call get_param(PF, mdl, "ICE_FLOAT_FRAC_VARNAME", floatfr_varname, &
                 "The name of the ice float fraction (grounding fraction) variable in ICE_VELOCITY_FILE.", &
                 default="float_frac")
  if (.not. skip_bed_local) then
    call get_param(PF, mdl, "BED_TOPO_FILE", bed_topo_file, &
                   "The file from which the bed elevation is read.", &
                   default="ice_shelf_vel.nc")
    call get_param(PF, mdl, "BED_TOPO_VARNAME", bed_varname, &
                   "The name of the bed elevation variable in ICE_INPUT_FILE.", &
                   default="depth")
  endif
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_ice_shelf_velocity_from_file: Unable to open "//trim(filename))

  call MOM_read_data(filename, trim(ushelf_varname), u_shelf, G%Domain, position=CORNER, scale=US%m_s_to_L_T)
  call MOM_read_data(filename, trim(vshelf_varname), v_shelf, G%Domain, position=CORNER, scale=US%m_s_to_L_T)
  call MOM_read_data(filename, trim(floatfr_varname), float_cond, G%Domain, scale=1.)

  if (.not. skip_bed_local) then
    filename = trim(inputdir)//trim(bed_topo_file)
    call MOM_read_data(filename, trim(bed_varname), bed_elev, G%Domain, scale=US%m_to_Z)
  endif


end subroutine initialize_ice_flow_from_file

!> Corner weights w(a,b) = int N(a,b) a(eta) d(xi) dxi deta, with d(xi) = dyW*(1-xi) + dyE*xi
!! and a(eta) = dxS*(1-eta) + dxN*eta. Each is a product of 1D integrals such as dyW/3 + dyE/6.
!! They sum to the cell area.
pure subroutine corner_cell_weights(dxS, dxN, dyW, dyE, w_cell)
  real, intent(in) :: dxS  !< South face length [L ~> m]
  real, intent(in) :: dxN  !< North face length [L ~> m]
  real, intent(in) :: dyW  !< West face length [L ~> m]
  real, intent(in) :: dyE  !< East face length [L ~> m]
  real, dimension(2,2), intent(out) :: w_cell !< Per-corner integration weights [L2 ~> m2]

  w_cell(1,1) = (dyW/3.0 + dyE/6.0) * (dxS/3.0 + dxN/6.0)
  w_cell(2,1) = (dyW/6.0 + dyE/3.0) * (dxS/3.0 + dxN/6.0)
  w_cell(1,2) = (dyW/3.0 + dyE/6.0) * (dxS/6.0 + dxN/3.0)
  w_cell(2,2) = (dyW/6.0 + dyE/3.0) * (dxS/6.0 + dxN/3.0)
end subroutine corner_cell_weights

!> Area-weighted cell mean of a bilinear field from its corner values and corner_cell_weights.
!! Opposite corners are paired for rotation invariance.
pure real function nodal_cell_mean(h_cell, w_cell) result(Hbar)
  real, dimension(2,2), intent(in) :: h_cell !< Corner values [A ~> a]
  real, dimension(2,2), intent(in) :: w_cell !< Per-corner integration weights [L2 ~> m2]

  real :: area  ! Sum of the corner weights, i.e. the cell area [L2 ~> m2]

  area = ((w_cell(1,1) + w_cell(2,2)) + (w_cell(1,2) + w_cell(2,1)))
  if (area > 0.0) then
    Hbar = ( ((w_cell(1,1)*h_cell(1,1)) + (w_cell(2,2)*h_cell(2,2))) + &
             ((w_cell(1,2)*h_cell(1,2)) + (w_cell(2,1)*h_cell(2,1))) ) / area
  else
    Hbar = 0.0
  endif
end function nodal_cell_mean

!> Area-weighted cell mean of the bilinear interpolant of the B-grid corner values of cell (i,j),
!! with the face lengths chosen as in init_nodal_DG_metric.
function corner_cell_mean(G, i, j, reentrant_x, reentrant_y, c11, c21, c12, c22) result(cmean)
  type(ocean_grid_type), intent(in) :: G   !< The grid structure used by the ice shelf.
  integer,               intent(in) :: i   !< The i-index of the cell.
  integer,               intent(in) :: j   !< The j-index of the cell.
  logical,               intent(in) :: reentrant_x !< True if the domain is zonally reentrant
  logical,               intent(in) :: reentrant_y !< True if the domain is meridionally reentrant
  real,                  intent(in) :: c11 !< Corner value at (I-1,J-1) [Z ~> m].
  real,                  intent(in) :: c21 !< Corner value at (I,J-1) [Z ~> m].
  real,                  intent(in) :: c12 !< Corner value at (I-1,J) [Z ~> m].
  real,                  intent(in) :: c22 !< Corner value at (I,J) [Z ~> m].
  real :: cmean                            !< The area-weighted cell mean [Z ~> m].

  real :: dxS, dxN, dyW, dyE      ! Face lengths [L ~> m]
  real, dimension(2,2) :: w_cell  ! Per-corner integration weights [L2 ~> m2]
  real, dimension(2,2) :: c_cell  ! Corner values in the (a,b) layout [Z ~> m]

  if ((J-1 >= G%JsdB) .and. (reentrant_y .or. (j + G%jdg_offset > G%jsg))) then
    dxS = G%dxCv(i,J-1) ; dxN = G%dxCv(i,J)
  else
    dxS = G%dxCv(i,J)   ; dxN = G%dxCv(i,J)
  endif
  if ((I-1 >= G%IsdB) .and. (reentrant_x .or. (i + G%idg_offset > G%isg))) then
    dyW = G%dyCu(I-1,j) ; dyE = G%dyCu(I,j)
  else
    dyW = G%dyCu(I,j)   ; dyE = G%dyCu(I,j)
  endif

  call corner_cell_weights(dxS, dxN, dyW, dyE, w_cell)
  c_cell(1,1) = c11 ; c_cell(2,1) = c21 ; c_cell(1,2) = c12 ; c_cell(2,2) = c22
  cmean = nodal_cell_mean(c_cell, w_cell)

end function corner_cell_mean

!> Read the nodal bed from BED_TOPO_FILE and set bed_elev to its cell means (INIT_ICE_BED_NODAL).
subroutine initialize_bed_node_from_file(bed_node, bed_elev, G, US, PF)
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: bed_node !< The nodal bed elevation [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: bed_elev !< The cell-averaged bed elevation [Z ~> m].
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type),  intent(in)    :: PF !< A structure to parse for run-time parameters

  character(len=200) :: filename, inputdir, nodal_bed_file
  character(len=200) :: bed_node_varname
  character(len=40)  :: mdl = "initialize_bed_node_from_file"
  logical :: reentrant_x, reentrant_y ! True if the domain wraps in x or y
  integer :: i, j

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "BED_TOPO_FILE", nodal_bed_file, &
                 "The file from which the nodal (B-grid corner) bed elevation is read "//&
                 "when INIT_ICE_BED_NODAL=True.", &
                 default="ice_shelf_vel.nc")
  call get_param(PF, mdl, "REENTRANT_X", reentrant_x, &
                 " If true, the domain is zonally reentrant.", &
                 default=.false., do_not_log=.true.)
  call get_param(PF, mdl, "REENTRANT_Y", reentrant_y, &
                 " If true, the domain is meridionally reentrant.", &
                 default=.false., do_not_log=.true.)
  call get_param(PF, mdl, "BED_TOPO_VARNAME", bed_node_varname, &
                 "The name of the nodal bed elevation variable in BED_TOPO_FILE "//&
                 "when INIT_ICE_BED_NODAL=True.", &
                 default="depth_n")

  filename = trim(inputdir)//trim(nodal_bed_file)
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_bed_node_from_file: Unable to open "//trim(filename))

  call MOM_read_data(filename, trim(bed_node_varname), bed_node, G%Domain, &
                     position=CORNER, scale=US%m_to_Z)
  call pass_var(bed_node, G%domain, position=CORNER)

  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    bed_elev(i,j) = corner_cell_mean(G, i, j, reentrant_x, reentrant_y, &
                                     bed_node(I-1,J-1), bed_node(I,J-1), &
                                     bed_node(I-1,J), bed_node(I,J))
  enddo ; enddo
  call pass_var(bed_elev, G%domain)

end subroutine initialize_bed_node_from_file

!> Initialize ice shelf b.c.s from file
subroutine initialize_ice_shelf_boundary_from_file(u_face_mask_bdry, v_face_mask_bdry, &
                u_bdry_val, v_bdry_val, umask, vmask, h_bdry_val, &
                hmask,  h_shelf, G, US, PF )

  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_face_mask_bdry !< A boundary-type mask at B-grid u faces [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_face_mask_bdry !< A boundary-type mask at B-grid v faces [nondim]
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: u_bdry_val !< The zonal ice shelf velocity at open
                                                      !! boundary vertices [L T-1 ~> m s-1].
  real, dimension(SZIB_(G),SZJB_(G)), &
                         intent(inout) :: v_bdry_val !< The meridional ice shelf velocity at open
                                                      !! boundary vertices [L T-1 ~> m s-1].
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: umask !< A mask for ice shelf velocity [nondim]
  real, dimension(SZDIB_(G),SZDJB_(G)), &
                         intent(inout) :: vmask !< A mask for ice shelf velocity [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: h_bdry_val !< The ice shelf thickness at open boundaries [Z ~> m]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: hmask !< A mask indicating which tracer points are
                                             !! partly or fully covered by an ice-shelf [nondim]
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(in) :: h_shelf !< Ice-shelf thickness [Z ~> m]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  character(len=200) :: filename, bc_file, inputdir, icethick_file ! Strings for file/path
  character(len=200) :: ufcmskbdry_varname, vfcmskbdry_varname, &
                        ubdryv_varname, vbdryv_varname, umask_varname, vmask_varname, &
                        hmsk_varname  ! Variable name in file
  character(len=40)  :: mdl = "initialize_ice_shelf_boundary_from_file" ! This subroutine's name.

  integer :: i, j, isc, jsc, iec, jec

  h_bdry_val(:,:) = 0.

  call MOM_mesg("  MOM_ice_shelf_init_profile.F90, initialize_b_c_s_from_file: reading b.c.s")

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "ICE_SHELF_BC_FILE", bc_file, &
                 "The file from which the boundary conditions are read.", &
                 default="ice_shelf_bc.nc")
  call get_param(PF, mdl, "ICE_THICKNESS_FILE", icethick_file, &
                 "The file from which the ice-shelf thickness is read.", &
                 default="ice_shelf_thick.nc")
  call get_param(PF, mdl, "ICE_THICKNESS_MASK_VARNAME", hmsk_varname, &
                 "The name of the icethickness mask variable in ICE_THICKNESS_FILE.", &
                 default="h_mask")

  filename = trim(inputdir)//trim(bc_file)
  call log_param(PF, mdl, "INPUTDIR/ICE_SHELF_BC_FILE", filename)
  call get_param(PF, mdl, "ICE_UBDRYMSK_VARNAME", ufcmskbdry_varname, &
                 "The name of the ice-shelf ubdrymask variable in ICE_SHELF_BC_FILE.", &
                 default="ufacemask")
  call get_param(PF, mdl, "ICE_VBDRYMSK_VARNAME", vfcmskbdry_varname, &
                 "The name of the ice-shelf vbdrymask variable in ICE_SHELF_BC_FILE.", &
                 default="vfacemask")
  call get_param(PF, mdl, "ICE_UMASK_VARNAME", umask_varname, &
                 "The name of the ice-shelf ubdrymask variable in ICE_SHELF_BC_FILE.", &
                 default="umask")
  call get_param(PF, mdl, "ICE_VMASK_VARNAME", vmask_varname, &
                 "The name of the ice-shelf vbdrymask variable in ICE_SHELF_BC_FILE.", &
                 default="vmask")
  call get_param(PF, mdl, "ICE_UBDRYVAL_VARNAME", ubdryv_varname, &
                 "The name of the ice-shelf ice_shelf ubdry variable in ICE_SHELF_BC_FILE.", &
                 default="ubdry_val")
  call get_param(PF, mdl, "ICE_VBDRYVAL_VARNAME", vbdryv_varname, &
                 "The name of the ice-shelf ice_shelf vbdry variable in ICE_SHELF_BC_FILE.", &
                 default="vbdry_val")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_ice_shelf_velocity_from_file: Unable to open "//trim(filename))


  call MOM_read_data(filename, trim(ufcmskbdry_varname), u_face_mask_bdry, G%Domain, position=CORNER, &
                     scale=1.)
  call MOM_read_data(filename, trim(vfcmskbdry_varname), v_face_mask_bdry, G%Domain, position=CORNER, &
                     scale=1.)
  call MOM_read_data(filename, trim(ubdryv_varname), u_bdry_val, G%Domain, position=CORNER, scale=US%m_s_to_L_T)
  call MOM_read_data(filename, trim(vbdryv_varname), v_bdry_val, G%Domain, position=CORNER, scale=US%m_s_to_L_T)
  call MOM_read_data(filename, trim(umask_varname), umask, G%Domain, position=CORNER, scale=1.)
  call MOM_read_data(filename, trim(vmask_varname), vmask, G%Domain, position=CORNER, scale=1.)
  filename = trim(inputdir)//trim(icethick_file)

  call MOM_read_data(filename,trim(hmsk_varname), hmask, G%Domain, scale=1.)
  isc = G%isc ; jsc = G%jsc ; iec = G%iec ; jec = G%jec

  do j=jsc,jec
    do i=isc,iec
      if (hmask(i,j) == 3.) then
        h_bdry_val(i,j) = h_shelf(i,j)
      endif
    enddo
  enddo

end subroutine initialize_ice_shelf_boundary_from_file

!> Initialize ice basal friction
subroutine initialize_ice_C_basal_friction(C_basal_friction, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: C_basal_friction !< Ice-stream basal friction
                                             !! in units of [R L Z T-2 (s m-1)^n_basal_fric ~> Pa (s m-1)^n_basal_fric]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

!  integer :: i, j
  real :: C_friction  ! Constant ice-stream basal friction in units of
                      ! [R L Z T-2 (s m-1)^n_basal_fric ~> Pa (s m-1)^n_basal_fric]
  character(len=40)  :: mdl = "initialize_ice_basal_friction" ! This subroutine's name.
  character(len=200) :: config
  character(len=200) :: varname
  character(len=200) :: inputdir, filename, C_friction_file

  call get_param(PF, mdl, "ICE_BASAL_FRICTION_CONFIG", config, &
                 "This specifies how the initial basal friction profile is specified. "//&
                 "Valid values are: CONSTANT and FILE.", &
                 fail_if_missing=.true.)

  if (trim(config)=="CONSTANT") then
    call get_param(PF, mdl, "BASAL_FRICTION_COEFF", C_friction, &
                 "Coefficient in sliding law.", units="Pa (s m-1)^(n_basal_fric)", default=5.e10, scale=US%Pa_to_RLZ_T2)

    C_basal_friction(:,:) = C_friction
  elseif (trim(config)=="FILE") then
    call MOM_mesg("  MOM_ice_shelf.F90, initialize_ice_shelf: reading friction coefficients")
    call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)

    call get_param(PF, mdl, "BASAL_FRICTION_FILE", C_friction_file, &
                "The file from which basal friction coefficients are read.", &
                default="ice_basal_friction.nc")
    filename = trim(inputdir)//trim(C_friction_file)
    call log_param(PF, mdl, "INPUTDIR/BASAL_FRICTION_FILE", filename)

    call get_param(PF, mdl, "BASAL_FRICTION_VARNAME", varname, &
                   "The variable to use in basal traction.", &
                   default="tau_b_beta")

    if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_ice_basal_friction_from_file: Unable to open "//trim(filename))

    call MOM_read_data(filename, trim(varname), C_basal_friction, G%Domain, scale=US%Pa_to_RLZ_T2)

  endif
end subroutine


!> Initialize ice-stiffness parameter
subroutine initialize_ice_AGlen(AGlen, ice_viscosity_compute, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: AGlen !< The ice-stiffness parameter A_Glen, often in [Pa-3 s-1]
  character(len=40) :: ice_viscosity_compute !< Specifies whether the ice viscosity is computed internally
                                             !! according to Glen's flow law; is constant (for debugging purposes)
                                             !! or using observed strain rates and read from a file
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  real :: A_Glen  ! Ice-stiffness parameter, often in [Pa-3 s-1]
  character(len=40)  :: mdl = "initialize_ice_stiffness" ! This subroutine's name.
  character(len=200) :: config
  character(len=200) :: varname
  character(len=200) :: inputdir, filename, AGlen_file

  call get_param(PF, mdl, "ICE_A_GLEN_CONFIG", config, &
                 "This specifies how the initial ice-stiffness parameter is specified. "//&
                 "Valid values are: CONSTANT and FILE.", &
                 fail_if_missing=.true.)

  if (trim(config)=="CONSTANT") then
    call get_param(PF, mdl, "A_GLEN", A_Glen, &
                   "Ice-stiffness parameter.", units="Pa-n_g s-1", default=2.261e-25)

    AGlen(:,:) = A_Glen

  elseif (trim(config)=="FILE") then
    call MOM_mesg("  MOM_ice_shelf.F90, initialize_ice_shelf: reading ice-stiffness parameter")
    call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)

    call get_param(PF, mdl, "ICE_STIFFNESS_FILE", AGlen_file, &
                 "The file from which the ice-stiffness is read.", &
                 default="ice_AGlen.nc")
    filename = trim(inputdir)//trim(AGlen_file)
    call log_param(PF, mdl, "INPUTDIR/ICE_STIFFNESS_FILE", filename)
    call get_param(PF, mdl, "A_GLEN_VARNAME", varname, &
                   "The variable to use as ice-stiffness.", &
                   default="A_GLEN")

    if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_ice_stiffness_from_file: Unable to open "//trim(filename))

    if (trim(ice_viscosity_compute) == "OBS") then
      ! AGlen is the ice viscosity [R L2 T-1 ~> Pa s] computed from obs and read from a file
      call MOM_read_data(filename, trim(varname), AGlen, G%Domain, scale=US%Pa_to_RL2_T2*US%s_to_T)
    else
      ! AGlen is the ice stiffness parameter [Pa-n_g s-1]
      call MOM_read_data(filename, trim(varname), AGlen, G%Domain)
    endif
  endif
end subroutine initialize_ice_AGlen

!> Initialize ice surface mass balance field that is held constant over time
subroutine initialize_ice_SMB(SMB, G, US, PF)
  type(ocean_grid_type), intent(in)    :: G    !< The ocean's grid structure
  real, dimension(SZDI_(G),SZDJ_(G)), &
                         intent(inout) :: SMB !< Ice surface mass balance parameter, often in [R Z T-1 ~> kg m-2 s-1]
  type(unit_scale_type), intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type), intent(in)    :: PF !< A structure to parse for run-time parameters

  real :: SMB_val  ! Constant ice surface mass balance parameter, often in [R Z T-1 ~> kg m-2 s-1]
  character(len=40)  :: mdl = "initialize_ice_SMB" ! This subroutine's name.
  character(len=200) :: config
  character(len=200) :: varname
  character(len=200) :: inputdir, filename, SMB_file

  call get_param(PF, mdl, "ICE_SMB_CONFIG", config, &
                 "This specifies how the initial ice surface mass balance parameter is specified. "//&
                 "Valid values are: CONSTANT and FILE.", &
                 default="CONSTANT")

  if (trim(config)=="CONSTANT") then
    call get_param(PF, mdl, "SMB", SMB_val, &
                 "Surface mass balance.", units="kg m-2 s-1", default=0.0, scale=US%kg_m2s_to_RZ_T)

    SMB(:,:) = SMB_val

  elseif (trim(config)=="FILE") then
    call MOM_mesg("  MOM_ice_shelf.F90, initialize_ice_shelf: reading SMB parameter")
    call get_param(PF, mdl, "INPUTDIR", inputdir, default=".")
    inputdir = slasher(inputdir)

    call get_param(PF, mdl, "ICE_SMB_FILE", SMB_file, &
                 "The file from which the ice surface mass balance is read.", &
                 default="ice_SMB.nc")
    filename = trim(inputdir)//trim(SMB_file)
    call log_param(PF, mdl, "INPUTDIR/ICE_SMB_FILE", filename)
    call get_param(PF, mdl, "ICE_SMB_VARNAME", varname, &
                   "The variable to use as surface mass balance.", &
                   default="SMB")

    if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_ice_SMV_from_file: Unable to open "//trim(filename))
    call MOM_read_data(filename,trim(varname), SMB, G%Domain, scale=US%kg_m2s_to_RZ_T)

  endif
end subroutine initialize_ice_SMB
end module MOM_ice_shelf_initialize
