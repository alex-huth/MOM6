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
public initialize_DG_thickness_slopes_from_file
public initialize_DG_thickness_from_node_file
public apply_DG1_inverse_mass
public project_corners_to_DG1_modal
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
subroutine initialize_ice_thickness(h_shelf, area_shelf_h, hmask, melt_mask, G, G_in, US, PF, rotate_index, turns)
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

  character(len=40)  :: mdl = "initialize_ice_thickness" ! This subroutine's name.
  character(len=200) :: config
  logical :: rotate = .false.
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
    allocate(tmp1_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=0.0)
    allocate(tmp2_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=0.0)
    allocate(tmp3_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=0.0)
    allocate(tmp4_2d(G_in%isd:G_in%ied,G_in%jsd:G_in%jed), source=1.0)
    select case ( trim(config) )
      case ("CHANNEL") ; call initialize_ice_thickness_channel (tmp1_2d, tmp2_2d, tmp3_2d, G_in, US, PF)
      case ("FILE") ; call initialize_ice_thickness_from_file (tmp1_2d, tmp2_2d, tmp3_2d, tmp4_2d, G_in, US, PF)
      case ("USER") ; call USER_init_ice_thickness (tmp1_2d, tmp2_2d, tmp3_2d, G_in, US, PF)
      case default  ; call MOM_error(FATAL,"MOM_initialize: Unrecognized ice profile setup "//trim(config))
    end select
    call rotate_array(tmp1_2d,turns, h_shelf)
    call rotate_array(tmp2_2d,turns, area_shelf_h)
    call rotate_array(tmp3_2d,turns, hmask)
    call rotate_array(tmp4_2d,turns, melt_mask)
    deallocate(tmp1_2d,tmp2_2d,tmp3_2d)
  else
    select case ( trim(config) )
      case ("CHANNEL") ; call initialize_ice_thickness_channel (h_shelf, area_shelf_h, hmask, G, US, PF)
      case ("FILE") ; call initialize_ice_thickness_from_file (h_shelf, area_shelf_h, hmask, melt_mask, G, US, PF)
      case ("USER") ; call USER_init_ice_thickness (h_shelf, area_shelf_h, hmask, G, US, PF)
      case default  ; call MOM_error(FATAL,"MOM_initialize: Unrecognized ice profile setup "//trim(config))
    end select
  endif

end subroutine initialize_ice_thickness

!> Initialize ice shelf thickness from file
subroutine initialize_ice_thickness_from_file(h_shelf, area_shelf_h, hmask, melt_mask, G, US, PF)
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

  !  This subroutine reads ice thickness and area from a file and puts it into
  !  h_shelf [Z ~> m] and area_shelf_h [L2 ~> m2] (and dimensionless) and updates hmask
  character(len=200) :: filename,thickness_file,inputdir ! Strings for file/path
  character(len=200) :: thickness_varname, area_varname, hmask_varname, melt_mask_varname  ! Variable name in file
  character(len=40)  :: mdl = "initialize_ice_thickness_from_file" ! This subroutine's name.
  integer :: i, j, isc, jsc, iec, jec
  logical :: hmask_set
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
                 "The name of the thickness variable in ICE_THICKNESS_FILE.", &
                 default="h_shelf")
  call get_param(PF, mdl, "ICE_AREA_VARNAME", area_varname, &
                 "The name of the area variable in ICE_THICKNESS_FILE.", &
                 default="area_shelf_h")
  hmask_varname="h_mask"
  call get_param(PF, mdl, "MELT_MASK_VARNAME", melt_mask_varname, &
                 "The name of the melt mask variable in ICE_THICKNESS_FILE.", &
                 default="melt_mask")
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_topography_from_file: Unable to open "//trim(filename))
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
    do j=jsc,jec ; do i=isc,iec
      hmask(i,j) = 0.0
      if (h_shelf(i,j) > 0.0) hmask(i,j) = 1.0
    enddo ; enddo
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
  logical, optional,     intent(in)    :: skip_bed !< If true, skip the BED_TOPO_FILE read.
                                                !! Used when bed_elev will be derived from a
                                                !! separately-read nodal bed field.

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

!> Read node-based bed elevation from NODAL_BED_FILE into CS%bed_node, then
!! derive the cell-centered CS%bed_elev by bilinear averaging the four
!! surrounding nodes. Replaces the reconstruct_bed_to_nodes Jacobi inversion
!! when an authoritative nodal bed product is available.
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
  real :: a0, a1, d0, d1   ! Per-cell metric scalars [L ~> m].
  real :: c1_unused, c2_unused ! Discarded modal slopes.
  integer :: i, j

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "BED_TOPO_FILE", nodal_bed_file, &
                 "The file from which the nodal (B-grid corner) bed elevation is read "//&
                 "when USE_NODAL_BED_FILE=True.", &
                 default="ice_shelf_vel.nc")
  call get_param(PF, mdl, "BED_TOPO_VARNAME", bed_node_varname, &
                 "The name of the nodal bed elevation variable in NODAL_BED_FILE.", &
                 default="depth_n")

  filename = trim(inputdir)//trim(nodal_bed_file)
  if (.not.file_exists(filename, G%Domain)) call MOM_error(FATAL, &
       " initialize_bed_node_from_file: Unable to open "//trim(filename))

  call MOM_read_data(filename, trim(bed_node_varname), bed_node, G%Domain, &
                     position=CORNER, scale=US%m_to_Z)
  call pass_var(bed_node, G%domain, position=CORNER)

  ! Derive cell-centered bed_elev as the area-weighted cell mean of the
  ! bilinear corner-node interpolant on the separable-Jacobian element.
  ! Uniform cells collapse to bed_elev = 0.25*sum(corners) bit-exactly.
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (J>1) then
      a0 = 0.5*(G%dxCv(i,J-1) + G%dxCv(i,J))
      a1 = G%dxCv(i,J) - G%dxCv(i,J-1)
    else
      a0 = G%dxCv(i,J) ; a1 = 0.0
    endif
    if (I>1) then
      d0 = 0.5*(G%dyCu(I-1,j) + G%dyCu(I,j))
      d1 = G%dyCu(I,j) - G%dyCu(I-1,j)
    else
      d0 = G%dyCu(I,j) ; d1 = 0.0
    endif
    call project_corners_to_DG1_modal(a0, a1, d0, d1, &
                                      bed_node(I-1,J-1), bed_node(I,J-1), &
                                      bed_node(I-1,J  ), bed_node(I,J  ), &
                                      bed_elev(i,j), c1_unused, c2_unused)
  enddo ; enddo
  call pass_var(bed_elev, G%domain)

end subroutine initialize_bed_node_from_file

!> Optionally read DG(1) cell-mean thickness slopes h_x, h_y from
!! ICE_THICKNESS_FILE. Returns slopes_set=.true. when both DG_HX_VARNAME and
!! DG_HY_VARNAME are present in the file; otherwise leaves h_x, h_y untouched
!! and returns slopes_set=.false. so the caller can fall back to its own
!! seeding (e.g. centred-FD of neighbour means).
subroutine initialize_DG_thickness_slopes_from_file(h_x, h_y, slopes_set, G, US, PF)
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_x !< Cell-mean DG(1) x-slope DOF [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_y !< Cell-mean DG(1) y-slope DOF [Z ~> m].
  logical,                intent(out)   :: slopes_set !< True if h_x and h_y were read from file.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type),  intent(in)    :: PF !< A structure to parse for run-time parameters

  character(len=200) :: filename, inputdir, thickness_file
  character(len=200) :: hx_varname, hy_varname
  character(len=40)  :: mdl = "initialize_DG_thickness_slopes_from_file"

  slopes_set = .false.

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "ICE_THICKNESS_FILE", thickness_file, &
                 default="ice_shelf_h.nc", do_not_log=.true.)
  call get_param(PF, mdl, "DG_HX_VARNAME", hx_varname, &
                 "The name of the DG(1) cell-mean x-slope variable in "//&
                 "ICE_THICKNESS_FILE. If both DG_HX_VARNAME and DG_HY_VARNAME "//&
                 "fields are present in the file, they are used as initial DG "//&
                 "thickness slopes (skipping the centred-FD seed and slope limiter, "//&
                 "mirroring restart behaviour).", &
                 default="h_x")
  call get_param(PF, mdl, "DG_HY_VARNAME", hy_varname, &
                 "The name of the DG(1) cell-mean y-slope variable in ICE_THICKNESS_FILE.", &
                 default="h_y")

  filename = trim(inputdir)//trim(thickness_file)
  if (.not.file_exists(filename, G%Domain)) return
  if (.not.field_exists(filename, trim(hx_varname), MOM_domain=G%Domain)) return
  if (.not.field_exists(filename, trim(hy_varname), MOM_domain=G%Domain)) return

  call MOM_read_data(filename, trim(hx_varname), h_x, G%Domain, scale=US%m_to_Z)
  call MOM_read_data(filename, trim(hy_varname), h_y, G%Domain, scale=US%m_to_Z)
  slopes_set = .true.

end subroutine initialize_DG_thickness_slopes_from_file

!> Optionally read nodal (B-grid corner) ice thickness from ICE_THICKNESS_FILE
!! and use it to derive both the cell-mean h_shelf and the DG(1) cell-mean
!! slopes h_x, h_y via a bilinear node fit. Returns used=.true. when the
!! nodal field is present and consumed; otherwise leaves h_shelf, h_x, h_y
!! untouched so the caller can fall back to the cell-mean read and centred-FD
!! slope seed. The nodal field is consumed locally and not retained.
subroutine initialize_DG_thickness_from_node_file(h_shelf, h_x, h_y, used, G, US, PF)
  type(ocean_grid_type),  intent(in)    :: G  !< The grid structure used by the ice shelf.
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_shelf !< Cell-mean ice thickness [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_x !< Cell-mean DG(1) x-slope DOF [Z ~> m].
  real, dimension(SZDI_(G),SZDJ_(G)), &
                          intent(inout) :: h_y !< Cell-mean DG(1) y-slope DOF [Z ~> m].
  logical,                intent(out)   :: used !< True if the nodal IC was applied.
  type(unit_scale_type),  intent(in)    :: US !< A structure containing unit conversion factors
  type(param_file_type),  intent(in)    :: PF !< A structure to parse for run-time parameters

  real, allocatable :: h_node(:,:) ! Nodal (corner) ice thickness [Z ~> m].
  character(len=200) :: filename, inputdir, thickness_file
  character(len=200) :: node_varname
  character(len=40)  :: mdl = "initialize_DG_thickness_from_node_file"
  real :: a0, a1, d0, d1   ! Per-cell metric scalars [L ~> m].
  integer :: i, j

  used = .false.

  call get_param(PF, mdl, "INPUTDIR", inputdir, default=".", do_not_log=.true.)
  inputdir = slasher(inputdir)
  call get_param(PF, mdl, "ICE_THICKNESS_FILE", thickness_file, &
                 default="ice_shelf_h.nc", do_not_log=.true.)
  call get_param(PF, mdl, "ICE_THICKNESS_NODE_VARNAME", node_varname, &
                 "The name of the nodal (B-grid corner) ice thickness variable in "//&
                 "ICE_THICKNESS_FILE. If present and USE_DG_THICKNESS=True, the field "//&
                 "is read and used to derive both the cell-mean h_shelf and the DG(1) "//&
                 "thickness slopes h_x, h_y via a bilinear node fit, skipping the "//&
                 "centred-FD slope seed and slope limiter. The nodal field is not "//&
                 "retained after initialization.", &
                 default="h_shelf_node")

  filename = trim(inputdir)//trim(thickness_file)
  if (.not.file_exists(filename, G%Domain)) return
  if (.not.field_exists(filename, trim(node_varname), MOM_domain=G%Domain)) return

  allocate(h_node(G%IsdB:G%IedB, G%JsdB:G%JedB)) ; h_node(:,:) = 0.0
  call MOM_read_data(filename, trim(node_varname), h_node, G%Domain, &
                     position=CORNER, scale=US%m_to_Z)
  call pass_var(h_node, G%domain, position=CORNER)

  ! Exact projection of the bilinear node interpolant onto the DG(1) monomial
  ! basis {1, xi, eta} with weight a(eta)*d(xi). Uniform cells reduce to the
  ! 0.25/0.5/0.5 formulas bit-exactly.
  do j=G%jsc,G%jec ; do i=G%isc,G%iec
    if (J>1) then
      a0 = 0.5*(G%dxCv(i,J-1) + G%dxCv(i,J))
      a1 = G%dxCv(i,J) - G%dxCv(i,J-1)
    else
      a0 = G%dxCv(i,J) ; a1 = 0.0
    endif
    if (I>1) then
      d0 = 0.5*(G%dyCu(I-1,j) + G%dyCu(I,j))
      d1 = G%dyCu(I,j) - G%dyCu(I-1,j)
    else
      d0 = G%dyCu(I,j) ; d1 = 0.0
    endif
    call project_corners_to_DG1_modal(a0, a1, d0, d1, &
                                      h_node(I-1,J-1), h_node(I,J-1), &
                                      h_node(I-1,J  ), h_node(I,J  ), &
                                      h_shelf(i,j), h_x(i,j), h_y(i,j))
  enddo ; enddo

  deallocate(h_node)
  used = .true.

end subroutine initialize_DG_thickness_from_node_file

!> Apply M_monomial^{-1} to a 3-vector of DG(1) right-hand-sides via the
!! tensor-product Gram-Schmidt shifted-orthogonal basis. M is the DG(1) mass
!! matrix on a separable-Jacobian element with weight a(eta)*d(xi), where
!! a(eta) = a0 + a1*eta and d(xi) = d0 + d1*xi on xi,eta in [-1/2, 1/2].
!! Input rhs0,rhs1,rhs2 are the unscaled volume-integral RHS components in
!! the monomial basis {1, xi, eta}. Output c0,c1,c2 are the corresponding
!! monomial coefficients of the DG polynomial. On uniform cells (a1=d1=0)
!! this collapses bit-exactly to c0 = rhs0/(a0*d0), c1 = 12*rhs1/(a0*d0),
!! c2 = 12*rhs2/(a0*d0).
pure subroutine apply_DG1_inverse_mass(a0, a1, d0, d1, rhs0, rhs1, rhs2, c0, c1, c2)
  real, intent(in)  :: a0   !< Cell-mean x metric (dxCv_S + dxCv_N)/2 [L ~> m]
  real, intent(in)  :: a1   !< Cell x-metric anisotropy dxCv_N - dxCv_S [L ~> m]
  real, intent(in)  :: d0   !< Cell-mean y metric (dyCu_W + dyCu_E)/2 [L ~> m]
  real, intent(in)  :: d1   !< Cell y-metric anisotropy dyCu_E - dyCu_W [L ~> m]
  real, intent(in)  :: rhs0 !< Volume-integral RHS for the {1} basis [Z L2 T-1 ~> m3 s-1]
  real, intent(in)  :: rhs1 !< Volume-integral RHS for the {xi} basis [Z L2 T-1 ~> m3 s-1]
  real, intent(in)  :: rhs2 !< Volume-integral RHS for the {eta} basis [Z L2 T-1 ~> m3 s-1]
  real, intent(out) :: c0   !< Monomial coefficient of 1 [Z T-1 ~> m s-1]
  real, intent(out) :: c1   !< Monomial coefficient of xi [Z T-1 ~> m s-1]
  real, intent(out) :: c2   !< Monomial coefficient of eta [Z T-1 ~> m s-1]

  real :: alpha_1, alpha_2 ! Shift parameters of the orthogonal basis [nondim].
  real :: Mtilde_00, Mtilde_11, Mtilde_22 ! Diagonal entries of the shifted mass matrix.
  real :: rhs0p, rhs1p, rhs2p ! Shifted-basis RHS components.
  real :: c0p, c1p, c2p ! Shifted-basis coefficients.

  alpha_1 = d1 / (12.0 * d0)
  alpha_2 = a1 / (12.0 * a0)

  Mtilde_00 = a0 * d0
  Mtilde_11 = a0 * ((12.0*d0*d0) - (d1*d1)) / (144.0 * d0)
  Mtilde_22 = d0 * ((12.0*a0*a0) - (a1*a1)) / (144.0 * a0)

  ! Forward affine transform to shifted basis.
  rhs0p = rhs0
  rhs1p = rhs1 - (alpha_1 * rhs0)
  rhs2p = rhs2 - (alpha_2 * rhs0)

  ! Diagonal scale.
  c0p = rhs0p / Mtilde_00
  c1p = rhs1p / Mtilde_11
  c2p = rhs2p / Mtilde_22

  ! Back transform to monomial coefficients.
  c0 = c0p - ((alpha_1 * c1p) + (alpha_2 * c2p))
  c1 = c1p
  c2 = c2p

end subroutine apply_DG1_inverse_mass

!> Exact projection of bilinear corner-node values onto DG(1) monomial
!! coefficients on a separable-Jacobian element. Corner values are at the four
!! reference corners (xi, eta) = (-1/2, -1/2), (+1/2, -1/2), (-1/2, +1/2),
!! (+1/2, +1/2) for (SW, SE, NW, NE). On uniform cells (a1=d1=0) this collapses
!! bit-exactly to c0 = 0.25*sum(corners), c1 = 0.5*((SE+NE)-(SW+NW)),
!! c2 = 0.5*((NW+NE)-(SW+SE)).
pure subroutine project_corners_to_DG1_modal(a0, a1, d0, d1, SW, SE, NW, NE, c0, c1, c2)
  real, intent(in)  :: a0  !< Cell-mean x metric [L ~> m].
  real, intent(in)  :: a1  !< Cell x-metric anisotropy [L ~> m].
  real, intent(in)  :: d0  !< Cell-mean y metric [L ~> m].
  real, intent(in)  :: d1  !< Cell y-metric anisotropy [L ~> m].
  real, intent(in)  :: SW  !< Corner value at (xi,eta)=(-1/2,-1/2) [<value units>].
  real, intent(in)  :: SE  !< Corner value at (xi,eta)=(+1/2,-1/2) [<value units>].
  real, intent(in)  :: NW  !< Corner value at (xi,eta)=(-1/2,+1/2) [<value units>].
  real, intent(in)  :: NE  !< Corner value at (xi,eta)=(+1/2,+1/2) [<value units>].
  real, intent(out) :: c0  !< Monomial coefficient of 1.
  real, intent(out) :: c1  !< Monomial coefficient of xi.
  real, intent(out) :: c2  !< Monomial coefficient of eta.

  real :: Pxm, Pxp, Qxm, Qxp ! 1D xi-integrals of d(xi) and xi*d(xi).
  real :: Pym, Pyp, Qym, Qyp ! 1D eta-integrals of a(eta) and eta*a(eta).
  real :: rhs0, rhs1, rhs2

  Pxm = (0.5*d0) - (d1/12.0)
  Pxp = (0.5*d0) + (d1/12.0)
  Qxm = (-d0/12.0) + (d1/24.0)
  Qxp = (d0/12.0) + (d1/24.0)

  Pym = (0.5*a0) - (a1/12.0)
  Pyp = (0.5*a0) + (a1/12.0)
  Qym = (-a0/12.0) + (a1/24.0)
  Qyp = (a0/12.0) + (a1/24.0)

  ! Diagonal + off-diagonal pairing in each rhs so the projection is bit-exact
  ! under a 90 deg grid rotation (which permutes corners SW<->SE<->NE<->NW).
  rhs0 = ((SW*(Pxm*Pym)) + (NE*(Pxp*Pyp))) + ((SE*(Pxp*Pym)) + (NW*(Pxm*Pyp)))
  rhs1 = ((SW*(Qxm*Pym)) + (NE*(Qxp*Pyp))) + ((SE*(Qxp*Pym)) + (NW*(Qxm*Pyp)))
  rhs2 = ((SW*(Pxm*Qym)) + (NE*(Pxp*Qyp))) + ((SE*(Pxp*Qym)) + (NW*(Pxm*Qyp)))

  call apply_DG1_inverse_mass(a0, a1, d0, d1, rhs0, rhs1, rhs2, c0, c1, c2)

end subroutine project_corners_to_DG1_modal

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
