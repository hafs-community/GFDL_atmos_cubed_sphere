!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANYWARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************
#ifdef OVERLOAD_R4
#define _GET_VAR1 get_var1_real
#else
#define _GET_VAR1 get_var1_double
#endif

!>@brief The module fv_nwp_nudge contains routines for nudging
!! to input analyses.
!>note This module is currently not supported in UFS

module fv_nwp_nudge_mod

! <table>
! <tr>
!     <th>Module Name</th>
!     <th>Functions Included</th>
!   </tr>
!   <tr>
!     <td>constants_mod</td>
!     <td>pi=>pi_8, grav, rdgas, cp_air, kappa, cnst_radius =>radius</td>
!   </tr>
!   <tr>
!     <td>external_sst_mod</td>
!     <td>i_sst, j_sst, sst_ncep, sst_anom, forecast_mode</td>
!   </tr>
!   <tr>
!     <td>diag_manager_mod</td>
!     <td>register_diag_field, send_data</td>
!   </tr>
!   <tr>
!     <td>fms_mod</td>
!     <td>write_version_number, open_namelist_file, check_nml_error,
!         file_exist, close_file</td>
!   </tr>
!   <tr>
!     <td>fv_arrays_mod</td>
!     <td>fv_grid_type, fv_grid_bounds_type, fv_nest_type, R_GRID</td>
!   </tr>
!   <tr>
!     <td>fv_diagnostics_mod</td>
!     <td>prt_maxmin, fv_time</td>
!   </tr>
!   <tr>
!     <td>fv_grid_utils_mod</td>
!     <td>great_circle_dist, intp_great_circle,
!         latlon2xyz, vect_cross, normalize_vect</td>
!   </tr>
!   <tr>
!     <td>fv_mp_mod</td>
!     <td>nmp_reduce_sum, mp_reduce_min, mp_reduce_max, is_master</td>
!   </tr>
!   <tr>
!     <td>fv_mapz_mod</td>
!     <td>mappm</td>
!   </tr>
!   <tr>
!     <td>fv_timing_mod</td>
!     <td>timing_on, timing_off</td>
!   </tr>
!   <tr>
!     <td>mpp_mod</td>
!     <td>mpp_error, FATAL, stdlog, get_unit, mpp_pe</td>
!   </tr>
!   <tr>
!     <td>sim_nc_mod</td>
!     <td>open_ncfile, close_ncfile, get_ncdim1, get_var1_double,
!         get_var3_r4, get_var1_real</td>
!   </tr>
!   <tr>
!     <td>mpp_domains_mod</td>
!     <td>mpp_update_domains, domain2d</td>
!   </tr>
!   <tr>
!     <td>time_manager_mod</td>
!     <td>time_type,  get_time, get_date</td>
!   </tr>
!   <tr>
!     <td>tp_core_mod</td>
!     <td>copy_corners</td>
!   </tr>
! </table>

 use external_sst_mod,  only: i_sst, j_sst, sst_ncep, sst_anom, forecast_mode
 use diag_manager_mod,  only: register_diag_field, send_data
#ifdef OVERLOAD_R4
 use constantsR4_mod,     only: pi=>pi_8, grav, rdgas, cp_air, kappa, cnst_radius =>radius
#else
 use constants_mod,     only: pi=>pi_8, grav, rdgas, cp_air, kappa, cnst_radius =>radius
#endif
 use fms_mod,           only: write_version_number, check_nml_error
 use fms2_io_mod,       only: file_exists
 use mpp_mod,           only: mpp_error, FATAL, stdlog, get_unit, mpp_pe, input_nml_file
 use mpp_domains_mod,   only: mpp_update_domains, domain2d
 use time_manager_mod,  only: time_type,  get_time, get_date
 use platform_mod,      only: r4_kind, r8_kind

 use fv_grid_utils_mod, only: great_circle_dist, intp_great_circle
 use fv_grid_utils_mod, only: latlon2xyz, vect_cross, normalize_vect
 use fv_diagnostics_mod,only: prt_maxmin, fv_time
 use tp_core_mod,       only: copy_corners
 use fv_mapz_mod,       only: mappm
 use fv_mp_mod,         only: mp_reduce_sum, mp_reduce_min, mp_reduce_max, is_master
 use fv_timing_mod,     only: timing_on, timing_off

 use sim_nc_mod,        only: open_ncfile, close_ncfile, get_ncdim1, get_var1_double, &
                              get_var3_r4, get_var1_real
 use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type, fv_nest_type, R_GRID
#ifdef MULTI_GASES
 use multi_gases_mod,  only:  virq, virqd, vicpqd, num_gas
#endif
 use mpp_domains_mod,   only: mpp_get_boundary, &
                              mpp_get_ntile_count, mpp_get_pelist, &
                              mpp_get_compute_domains, mpp_global_field, &
                              mpp_get_data_domain, mpp_get_compute_domain, &
                              mpp_get_global_domain

 implicit none
 private

 real(kind=R_GRID), parameter :: radius = cnst_radius

! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

 logical :: do_adiabatic_init

 public fv_nwp_nudge, fv_nwp_nudge_init, fv_nwp_nudge_end, breed_slp_inline, T_is_Tv
 public do_adiabatic_init
 public nwp_nudge_int
 integer im     !< Data x-dimension
 integer jm     !< Data y-dimension
 integer km     !< Data z-dimension
 real, allocatable:: ak0(:), bk0(:)
 real, allocatable:: lat(:), lon(:)

 logical :: module_is_initialized = .false.
 logical :: master
 logical :: no_obs = .false.
 real :: deg2rad, rad2deg
 real :: time_nudge = 0.
 integer :: time_interval = 6*3600   !< dataset time interval (seconds)
! ---> h1g, enhance the max. analysis data files, 2012-10-22
! integer, parameter :: nfile_max = 125
 integer, parameter :: nfile_max = 29280  !< maximum: 20-year analysis data, 4*366*20=29280
! <--- h1g,  2012-10-22
 integer :: nfile

 integer :: k_breed = 0
 integer :: k_trop = 0
 real    :: p_trop = 950.E2
 real    :: dps_min = 50.      !< maximum PS increment (Pa; each step) due to inline breeding
 real    :: del2_cd = 0.16

 real,    allocatable:: s2c(:,:,:)
 integer, allocatable:: id1(:,:), id2(:,:), jdc(:,:)
 real, allocatable :: ps_dat(:,:,:)
 real(KIND=4), allocatable, dimension(:,:,:,:):: u_dat, v_dat, t_dat, q_dat, dp_dat
 real(KIND=4), allocatable, dimension(:,:,:):: gz3  !< work array
 real, allocatable:: gz0(:,:)

! Namelist variables:
! ---> h1g, add the list of input NCEP analysis data files, 2012-10-22
 character(len=128):: input_fname_list =""       !< a file lists the input NCEP analysis data
 character(len=128):: analysis_file_first =""    !< the first NCEP analysis file to be used for nudging,
                                                 !! by default, the first file in the "input_fname_list"
 character(len=128):: analysis_file_last=""      !< the last NCEP analysis file to be used for nudging
                                                 !! by default, the last file in the "input_fname_list"

 real   :: P_relax = 30.E2                       !< from P_relax upwards, nudging is reduced linearly
                                                 !! proportional to pfull/P_relax
 real   :: P_Lower = 99999.E2                   !< from P_Lower downwards, nudging is reduced gradually

 real   :: P_norelax = 0.0                       !< from P_norelax upwards, no nudging
! <--- h1g, 2012-10-22

 character(len=128):: file_names(nfile_max)
 character(len=128):: track_file_name
 integer :: nfile_total = 0       !< =5 for 1-day (if datasets are 6-hr apart)
 real    :: p_wvp = 100.E2        !< cutoff level for specific humidity nudging
 integer :: kord_data = 8

 integer :: nwp_nudge_int  = 21600   ! 6 hours

 real    :: mask_fac = 0.25       !< [0,1]  0: no mask;  1: full strength
 integer :: nudge_option = 0    !< [0,1,2,3,4] 0: original global nudging; 1: regional grid nudging; 2: regional spectural nudging; 3. global grid nudging; 4. global spectral nudging

 logical :: T_is_Tv  = .false.
 logical :: use_pt_inc  = .false.
 logical :: use_high_top = .false.
 logical :: add_bg_wind = .true.
 logical :: conserve_mom = .true.
 logical :: conserve_hgt = .true.
 logical :: tc_mask = .false.
 logical :: strong_mask = .false.
 logical :: ibtrack = .true.
 logical :: nudge_debug = .false.
 logical :: rf_debug_dump = .false.   ! dump selected SN slices before/after RF
 integer :: rf_debug_k = 50          ! model level used by rf_debug_dump
 logical :: rf_cube_continuous = .true. ! T: continuous 6-tile RF; F: original independent-tile RF
 logical :: do_ps_bias  = .false.
 logical :: nudge_ps    = .false.
 logical :: nudge_q     = .false.
 logical :: nudge_winds = .true.
 logical :: nudge_virt  = .true.
 logical :: nudge_hght  = .true.
 logical :: time_varying = .true.
 logical :: time_varying_nwp = .false.
 logical :: print_end_breed = .true.
 logical :: print_end_nudge = .true.
 logical :: using_merra2 = .false.  ! Flag to allow avoidance of multiplicative factor if using MERRA2 data.
 logical :: climate_nudging = .false. ! Flag to allow for climate nudging.

 ! Optional direct pressure-level GRIB2 input for nudge_option > 0.
 !     .false. keeps the existing get_regional_analysis NetCDF path unchanged.
 !     .true.  reads a regular lat-lon GRIB2 file directly, horizontally
 !             interpolates to the FV3 A-grid, and vertically interpolates
 !             the pressure-level point data to the actual FV3 full levels.
 !     grib2_nplev and grib2_plev_pa may be overridden in fv_nwp_nudge_nml.
 !     Pressure levels are in Pa and must increase downward.
 integer, parameter :: grib2_nplev_default = 13
 integer, parameter :: grib2_nplev_max = 100
 real, parameter :: grib2_default_plev_pa(grib2_nplev_default) = (/ &
      5000.0, 10000.0, 15000.0, 20000.0, 25000.0, 30000.0, 40000.0, &
     50000.0, 60000.0, 70000.0, 85000.0, 92500.0,100000.0 /)
 logical :: grib2_input = .false.
 integer :: grib2_nplev = grib2_nplev_default
 real :: grib2_plev_pa(grib2_nplev_max)
 ! GRIB2 source-grid transform saved from the first field so every later
 ! variable/level is canonicalized in exactly the same way.
 integer, save :: grib2_raw_im = 0, grib2_raw_jm = 0
 integer, save :: grib2_lon_start = 1
 logical, save :: grib2_drop_duplicate_lon = .false.
 logical, save :: grib2_grid_signature_set = .false.
 integer, save :: grib2_grid_signature(19) = 0
 ! Within one GRIB2 file, lock all fields to the same reference/forecast time
 ! as the first selected field so a multi-time file cannot silently mix times.
 logical, save :: grib2_time_match_set = .false.
 integer, save :: grib2_time_jids(7) = -9999
 integer, save :: grib2_time_unit = -9999, grib2_forecast_time = -9999

!Nudging time-scales (seconds):
!note, however, the effective time-scale is 2X smaller (stronger) due
! to the use of the time-varying weighting factor
 real :: tau_ps     = 21600.       !< 1-day
 real :: tau_q      = 86400.       !< 1-day
 real :: tau_winds  = 21600.       !< 6-hr
 real :: tau_virt   = 43200.
 real :: tau_hght   = 43200.

 real :: q_min      = 1.E-8

 integer :: jbeg, jend
 integer :: nf_uv = 0
 integer :: nf_ps = 0
 integer :: nf_t  = 2
 integer :: nf_ht = 1
 integer :: L_uv  = 120 ! RF length_scale
 integer :: L_t   = 120
 ! Independent latitude tapers for temperature and winds.
 !     <= 0.0 : disabled (weight = 1 everywhere)
 !     >  0.0 : full nudging poleward of this latitude, with an
 !              exponential taper toward the equator.
 real :: lat_nudge_taper_t  = 0.0
 real :: lat_nudge_taper_uv = 0.0
 real, allocatable, save :: lat_weight_t(:,:), lat_weight_uv(:,:)
 logical, save :: lat_taper_initialized = .false.
 integer :: nf_skip_timestep = 1
 integer :: nf_skip_grid   = 0
 integer :: skip_counter= 0

! starting layer (top layer is sponge layer and is skipped)
 integer :: kstart = 2

! skip "kbot" layers
 integer :: kbot_winds = 0
 integer :: kbot_t     = 0
 integer :: kbot_q     = 0
 logical :: analysis_time

!-- Tropical cyclones  --------------------------------------------------------------------

! track dataset: 'INPUT/tropical_cyclones.txt'

  logical :: breed_srf_w = .false.
  real :: grid_size     = 28.E3
  real :: tau_vt_slp    = 1200.
  real :: tau_vt_wind   = 1200.
  real :: tau_vt_rad    = 4.0

  real :: pt_lim =  0.2
  real ::  slp_env = 101010.    !< storm environment pressure (Pa)
  real :: pre0_env = 100000.    !< critical storm environment pressure (Pa) for size computation
  real, parameter:: tm_max = 315.
!------------------
  real:: r_lo = 2.0
  real:: r_hi = 5.0    ! try 4.0?
!------------------
  real::  r_fac = 1.2
  real :: r_min = 200.E3
  real :: r_inc =  25.E3
  real, parameter:: del_r = 50.E3
  real:: elapsed_time = 0.0
  real:: nudged_time = 1.E12 !< seconds
                             ! usage example: set to 43200. to do inline vortex breeding
                             ! for only the first 12 hours
                             ! In addition, specify only 3 analysis files (12 hours)
  integer:: year_track_data
  integer, parameter:: max_storm = 140     !< max # of storms to process
  integer, parameter::  nobs_max = 125     !< Max # of observations per storm

  integer :: nstorms = 0
  integer :: nobs_tc(max_storm)
  integer :: min_nobs = 16
  real :: min_mslp = 1009.E2
  real(KIND=4)::     x_obs(nobs_max,max_storm)           !< longitude in degrees
  real(KIND=4)::     y_obs(nobs_max,max_storm)           !< latitude in degrees
  real(KIND=4)::  wind_obs(nobs_max,max_storm)           !< observed 10-m wind speed (m/s)
  real(KIND=4)::  mslp_obs(nobs_max,max_storm)           !< observed SLP in mb
  real(KIND=4)::  mslp_out(nobs_max,max_storm)           !< outer ring SLP in mb
  real(KIND=4)::   rad_out(nobs_max,max_storm)           !< outer ring radius in meters
  real(KIND=4)::   time_tc(nobs_max,max_storm)           !< start time of the track
!------------------------------------------------------------------------------------------
  integer :: id_ht_err

  !!! CLEANUP Module pointers
  integer :: is, ie, js, je
  integer :: isd, ied, jsd, jed

  ! cache (global)
  integer, save :: nudge_tile_num = -1
  logical, save :: rf_alpha_stored = .false.
  integer, save :: rf_alpha_mode_stored = 0
  integer, save :: rf_isg = 0, rf_ieg = -1, rf_jsg = 0, rf_jeg = -1
  integer, save :: rf_alpha_stride = 1
  integer, save :: rf_alpha_tile_stored = -1, rf_alpha_root_stored = -1
  integer, save :: rf_alpha_isc_stored = 0, rf_alpha_iec_stored = -1
  integer, save :: rf_alpha_jsc_stored = 0, rf_alpha_jec_stored = -1
  real,    save :: rf_alpha_Lx_m_stored = -1.0
  real,    save :: rf_alpha_Ly_m_stored = -1.0
  real,    save :: rf_alpha_dx_min_stored = 0.0, rf_alpha_dx_max_stored = 0.0
  real,    save :: rf_alpha_dy_min_stored = 0.0, rf_alpha_dy_max_stored = 0.0
  real,    save :: rf_alpha_dx_sum_stored = 0.0, rf_alpha_dy_sum_stored = 0.0
  real,    save :: rf_alpha_dx_sq_sum_stored = 0.0, rf_alpha_dy_sq_sum_stored = 0.0
  real, allocatable, save :: alpha_x_glob(:,:), alpha_y_glob(:,:)
  ! Centralized cubed-sphere RF storage.  For nudge_option=4 all tiles
  ! are gathered to PE 0.  The third dimension is the FV3 tile number.
  real, allocatable, save :: alpha_x_cube(:,:,:), alpha_y_cube(:,:,:)
  integer, save :: rf_cube_ntiles = 0
  integer, save :: rf_neighbor_tile(6,4) = 0
  integer, save :: rf_neighbor_edge(6,4) = 0
  logical, save :: rf_neighbor_reverse(6,4) = .false.
  ! Cached continuous-cube path geometry on the active RF grid.
  logical, save :: rf_cube_paths_built = .false.
  integer, save :: rf_path_nx = 0, rf_path_ny = 0, rf_path_ntiles = 0
  integer, save :: rf_path_nstate = 0, rf_path_nloop = 0, rf_path_maxloop = 0
  integer, allocatable, save :: rf_path_t(:), rf_path_i(:), rf_path_j(:), rf_path_dir(:)
  integer, allocatable, save :: rf_loop_start(:), rf_loop_len(:)
  integer, allocatable, save :: rf_hit_count(:,:,:)
  real, allocatable, save :: dt_stored(:,:,:), du_stored(:,:,:), dv_stored(:,:,:)
  integer, save :: nact
  integer, allocatable, save :: klist(:)

 namelist /fv_nwp_nudge_nml/T_is_Tv, nudge_ps, nudge_virt, nudge_hght, nudge_q, nudge_winds, &
                          do_ps_bias, tau_ps, tau_winds, tau_q, tau_virt, tau_hght,  kstart, kbot_winds, &
                          k_breed, k_trop, p_trop, dps_min, kord_data, tc_mask, nudge_option, nudge_debug, &
                          rf_debug_dump, rf_debug_k, rf_cube_continuous, nf_ps, nf_t, nf_ht,  &
                          nf_uv, nf_skip_timestep, nf_skip_grid, breed_srf_w, tau_vt_wind, tau_vt_slp, strong_mask, mask_fac, del2_cd,   &
                          L_uv, L_t, lat_nudge_taper_t, lat_nudge_taper_uv, kbot_t, kbot_q, p_wvp, time_varying, time_interval, use_pt_inc, pt_lim,  &
                          tau_vt_rad, r_lo, r_hi, use_high_top, add_bg_wind, conserve_mom, conserve_hgt,  &
                          min_nobs, min_mslp, nudged_time, r_fac, r_min, r_inc, ibtrack, track_file_name, file_names,         &
                          input_fname_list, analysis_file_first, analysis_file_last, P_relax, P_Lower, P_norelax, &
                          nwp_nudge_int, time_varying_nwp, using_merra2, climate_nudging, grib2_input, &
                          grib2_nplev, grib2_plev_pa

 contains

!>@brief Ths subroutine 'fv_nwp_nudge' computes and returns time tendencies for nudging to analysis.
!>@details This nudging is typically applied to fv_update_phys.
  subroutine fv_nwp_nudge ( Time, dt, npx, npy, npz, ps_dt, u_dt, v_dt, t_dt, q_dt, zvir, ptop, &
                            ak, bk, ts, ps, delp, ua, va, pt, nwat, q, phis, gridstruct, &
                            bd, domain)

  type(time_type), intent(in):: Time
  integer,         intent(IN):: npx, npy
  integer,         intent(in):: npz           !< vertical dimension
  integer,         intent(in):: nwat
  real,            intent(in):: dt
  real,            intent(in):: zvir, ptop
  type(domain2d), intent(INOUT), target :: domain
  type(fv_grid_bounds_type), intent(IN) :: bd
  real, intent(in   ), dimension(npz+1):: ak, bk
  real, intent(in   ), dimension(isd:ied,jsd:jed    ):: phis
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: pt, ua, va, delp
! pt as input is true tempeature
#ifdef MULTI_GASES
  real, intent(inout):: q(isd:ied,jsd:jed,npz,*)
#else
  real, intent(inout):: q(isd:ied,jsd:jed,npz,nwat)
#endif
  real, intent(inout), dimension(isd:ied,jsd:jed):: ps
! Accumulated tendencies
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: u_dt, v_dt
  real, intent(out), dimension(is:ie,js:je,npz):: t_dt
  real, intent(inout), dimension(is:ie,js:je,npz,1):: q_dt
  real, intent(out), dimension(is:ie,js:je):: ps_dt, ts

  type(fv_grid_type), intent(INOUT), target :: gridstruct
! local:
  real:: h2(is:ie,js:je)
  real:: m_err(is:ie,js:je)         !< height error at specified model interface level
  real:: slp_n(is:ie,js:je)         !< "Observed" SLP
  real:: slp_m(is:ie,js:je)         !< "model" SLP
  real::   mask(is:ie,js:je)
  real::     tv(is:ie,js:je)
  real:: peln(is:ie,npz+1)
  real:: pe2(is:ie, npz+1)
  real:: ptmp
  real, allocatable :: ps_obs(:,:)
  real, allocatable, dimension(:,:,:):: u_obs, v_obs, t_obs, q_obs
  real, allocatable, dimension(:,:,:):: du_obs, dv_obs
  real:: ps_fac(is:ie,js:je)
  integer :: seconds, days
  integer :: i,j,k, iq, kht
  real :: factor, rms, bias, co
  real :: factor_nwp
  real :: rdt, press(npz), profile(npz), prof_t(npz), prof_q(npz), du, dv, profile_Lower(npz)
  real :: alat_deg, lat_diag, wt_t_diag, wt_uv_diag
  integer :: ilat
  logical used

  real :: Lx_m, Ly_m

  real(kind=R_GRID), pointer, dimension(:,:,:) :: agrid
  real, pointer, dimension(:,:)   :: rarea, area

  real, pointer, dimension(:,:) :: sina_u, sina_v
  real, pointer, dimension(:,:,:) :: sin_sg
  real(kind=R_GRID), pointer, dimension(:,:,:) :: vlon, vlat

  real, pointer, dimension(:,:) :: dx, dy, rdxc, rdyc

  real(kind=R_GRID), pointer :: da_min

  logical, pointer :: bounded_domain, sw_corner, se_corner, nw_corner, ne_corner

  integer :: c0, c1, rate
  real    :: tcost
  logical :: do_nudge_now
  if ( .not. module_is_initialized ) then
        call mpp_error(FATAL,'==> Error from fv_nwp_nudge: module not initialized')
  endif
  if (nudge_option .ne. 0) call system_clock(count_rate=rate)
  agrid => gridstruct%agrid_64
   area => gridstruct%area
  rarea => gridstruct%rarea

  vlon   => gridstruct%vlon
  vlat   => gridstruct%vlat
  sina_u => gridstruct%sina_u
  sina_v => gridstruct%sina_v
  sin_sg => gridstruct%sin_sg

  dx     => gridstruct%dx
  dy     => gridstruct%dy
  rdxc   => gridstruct%rdxc
  rdyc   => gridstruct%rdyc

  da_min => gridstruct%da_min

  sw_corner => gridstruct%sw_corner
  se_corner => gridstruct%se_corner
  nw_corner => gridstruct%nw_corner
  ne_corner => gridstruct%ne_corner

  if ( no_obs ) then
#ifndef DYCORE_SOLO
       forecast_mode = .true.
#endif
       return
  endif

  call get_time (time, seconds, days)

  do j=js,je
     do i=is,ie
        mask(i,j) = 1.
     enddo
  enddo
  if ( tc_mask )  call get_tc_mask(time, mask, agrid)

! Latitude-dependent nudging weights for T and U/V.
  if (nudge_option .ne. 0 .and. .not. lat_taper_initialized) then
     if (.not.allocated(lat_weight_t)) allocate(lat_weight_t(is:ie,js:je))
     if (.not.allocated(lat_weight_uv)) allocate(lat_weight_uv(is:ie,js:je))
     lat_weight_t(:,:)  = 1.0
     lat_weight_uv(:,:) = 1.0

     if (lat_nudge_taper_t > 0.0 .or. lat_nudge_taper_uv > 0.0) then
!$OMP parallel do default(none) &
!$OMP shared(is,ie,js,je,agrid,lat_weight_t,lat_weight_uv,lat_nudge_taper_t,lat_nudge_taper_uv) &
!$OMP private(i,j,alat_deg)
        do j=js,je
           do i=is,ie
              alat_deg = abs(agrid(i,j,2)) * 180.0 / pi
              if (lat_nudge_taper_t > 0.0 .and. alat_deg < lat_nudge_taper_t) then
                 lat_weight_t(i,j) = exp(-10.0 * &
                      ((lat_nudge_taper_t-alat_deg)/lat_nudge_taper_t)**2)
              endif
              if (lat_nudge_taper_uv > 0.0 .and. alat_deg < lat_nudge_taper_uv) then
                 lat_weight_uv(i,j) = exp(-10.0 * &
                      ((lat_nudge_taper_uv-alat_deg)/lat_nudge_taper_uv)**2)
              endif
           enddo
        enddo
!$OMP end parallel do
     endif
     lat_taper_initialized = .true.
  endif

! The following profile is suitable only for nwp purposes; if the analysis has a good representation
! of the strat-meso-sphere the profile for upper layers should be changed.

  profile(:) = 1.
  profile_Lower(:) = 1.

!$OMP parallel do default(none) shared(npz,press,ak,bk)
  do k=1,npz
     press(k) = 0.5*(ak(k) + ak(k+1)) + 0.5*(bk(k)+bk(k+1))*1.E5
  end do
!$OMP parallel do default(none) shared(npz,press,P_relax,P_norelax,profile,profile_Lower,P_Lower,nudge_option)
  do k=1,npz
     if ( press(k) < P_relax ) then
       if ( nudge_option .ne. 0 ) then
          profile(k) = max(0.01, exp(-10. * ((P_relax - press(k)) / P_relax)**2))
       else
          profile(k) =  max(0.01, press(k)/P_relax)
       endif
     endif
     if ( nudge_option .ne. 0 .and. press(k) > P_Lower ) then
        profile_Lower(k) = max(0.01, exp(-10. * ((press(k)-P_Lower)/(press(npz)-P_Lower))**2))
     endif
  ! above P_norelax, no nudging. added by h1g
     if( press(k) < P_norelax ) profile(k) = 0.0
  enddo
  profile(1) = 0.

! Thermodynamics:
  prof_t(:) = 1.
  if (climate_nudging) then
!$OMP parallel do default(none) shared(npz,press,P_norelax,prof_t)
    do k=1,npz
      if ( press(k) < 10.E2 ) then
        prof_t(k) =  max(0.01, press(k)/10.E2)
      endif
      ! above P_norelax, no nudging.
      if( press(k) < P_norelax ) prof_t(k) = 0.0
    enddo
  else
!$OMP parallel do default(none) shared(npz,press,prof_t)
    do k=1,npz
      if ( press(k) < 10.E2 ) then
          prof_t(k) =  max(0.01, press(k)/10.E2)
      endif
    enddo
  endif
  prof_t(1) = 0.
  if ( nudge_option .ne. 0 ) prof_t=profile ! Force T profile to follow wind

! Water vapor:
  prof_q(:) = 1.
  if ( climate_nudging) then
!$OMP parallel do default(none) shared(npz,press,P_norelax,prof_q)
    do k=1,npz
      if ( press(k) < 200.E2 ) then
          prof_q(k) =  max(0., press(k)/200.E2)
      endif
      if( press(k) < P_norelax ) prof_q(k) = 0.0
    enddo
  else
!$OMP parallel do default(none) shared(npz,press,prof_q)
    do k=1,npz
      if ( press(k) < 300.E2 ) then
          prof_q(k) =  max(0., press(k)/300.E2)
      endif
    enddo
  endif
  prof_q(1) = 0.

! Height
  if ( k_trop == 0 ) then
       k_trop = 2
       do k=2,npz-1
          ptmp = ak(k+1) + bk(k+1)*1.E5
          if ( ptmp > p_trop ) then
               k_trop = k
               exit
          endif
       enddo
  endif

  if ( nudge_virt .and. nudge_hght ) then
       kht = k_trop
  else
       kht = npz-kbot_t
  endif

  if ( time_varying ) then
       factor = 1. + cos(real(mod(seconds,time_interval))/real(time_interval)*2.*pi)
       factor = max(1.e-5, factor)
  else
       factor = 1.
  endif

  if ( time_varying_nwp ) then
    if (mod(seconds, nwp_nudge_int) == 0) then
      factor_nwp = 1.0
    else
      factor_nwp = 0.0
    endif
  else
    factor_nwp = 1.0
  endif

  if ( do_adiabatic_init ) factor = 2.*factor

  allocate (ps_obs(is:ie,js:je) )
  allocate ( t_obs(is:ie,js:je,npz) )
  allocate ( q_obs(is:ie,js:je,npz) )

  if ( nudge_winds ) then
       allocate (u_obs(is:ie,js:je,npz) )
       allocate (v_obs(is:ie,js:je,npz) )
  else if ( nudge_option .ne. 0 ) then
       ! New input paths still pass explicit-shape u_obs/v_obs dummies.
       allocate (u_obs(is:ie,js:je,npz) )
       allocate (v_obs(is:ie,js:je,npz) )
       u_obs = 0.0
       v_obs = 0.0
  endif


  call get_obs(Time, dt, zvir, ak, bk, ps, ts, ps_obs, delp, pt, nwat, q, u_obs, v_obs, t_obs, q_obs,   &
               phis, ua, va, u_dt, v_dt, npx, npy, npz, factor, factor_nwp, mask, ptop, bd, gridstruct, domain, nudge_tile_num)
! *t_obs* is virtual temperature for normal nudging and true temperature for climate_nudging.

  if ( no_obs ) then
       deallocate (ps_obs)
       deallocate (t_obs)
       deallocate (q_obs)
       if ( nudge_winds .or. nudge_option .ne. 0 ) then
            deallocate (u_obs)
            deallocate (v_obs)
       endif
#ifndef DYCORE_SOLO
       forecast_mode = .true.
#endif
       return
   endif

   if ( grib2_input ) then
      ps_fac(:,:) = 1.0
   else
      do j=js,je
         do i=is,ie
            if ( abs(ps(i,j)-ps_obs(i,j)) > 2.e2 ) then
                 ps_fac(i,j) = 2.e2 / abs(ps(i,j)-ps_obs(i,j))
            else
                 ps_fac(i,j) = 1.
            endif
         enddo
      enddo
   endif

  if (nudge_option .ne. 0 .and. nf_skip_timestep.gt.1) then
    if (skip_counter.eq.0) then
     skip_counter=nf_skip_timestep
     do_nudge_now=.true.
     if (.not. allocated(dt_stored))allocate(dt_stored(is:ie,js:je,1:npz))
     if (.not. allocated(du_stored))allocate(du_stored(is:ie,js:je,1:npz))
     if (.not. allocated(dv_stored))allocate(dv_stored(is:ie,js:je,1:npz))
    else
     do_nudge_now=.false.
    end if
  else
     do_nudge_now=.true.
  end if
  if (nudge_option .ne. 0) then
    nact=npz
    if (.not. allocated(klist)) allocate(klist(npz))
    do k=1,npz
       klist(k)=k
    enddo
    if (do_nudge_now) then
      nact = 0
      do k=1,npz
        if (abs(profile_Lower(k)) > 1.e-2) then
          nact = nact + 1
          klist(nact) = k
        end if
      end do
    end if
  endif
  if( analysis_time .and. .not. grib2_input ) then
!-------------------------------------------
! Compute RMSE, bias, and correlation of SLP
! Pressure-level GRIB2 has no analysis ps/terrain, so this diagnostic is skipped.
!-------------------------------------------
      do j=js,je
         do i=is,ie
#ifdef MULTI_GASES
            tv(i,j) = pt(i,j,npz)*virq(q(i,j,npz,1:num_gas))
#else
            tv(i,j) = pt(i,j,npz)*(1.+zvir*q(i,j,npz,1))
#endif
         enddo
      enddo
      call compute_slp(is, ie, js, je, tv, ps(is:ie,js:je), phis(is:ie,js:je), slp_m)
      if ( nudge_option.eq.0) then
       call compute_slp(is, ie, js, je, t_obs(is,js,npz), ps_obs, gz0, slp_n)
      else
       gz0(is:ie,js:je)=phis(is:ie,js:je) !Use phis from the forecast?
       call compute_slp(is, ie, js, je, t_obs(:,:,npz), ps_obs, gz0, slp_n)
      endif

      if ( nudge_debug ) then
           if(master) write(*,*) 'kht=', kht
           call prt_maxmin('PS_o', ps_obs, is, ie, js, je, 0, 1, 0.01)
           ptmp = 0.01*g0_sum(ps_obs, is, ie, js, je, 1, .false., isd, ied, jsd, jed, area)
           if(master) write(*,*) 'Mean PS_o=', ptmp
           call prt_maxmin('SLP_m', slp_m, is, ie, js, je, 0, 1, 0.01)
           call prt_maxmin('SLP_o', slp_n, is, ie, js, je, 0, 1, 0.01)
      endif

!$OMP parallel do default(none) shared(is,ie,js,je,phis,gz0,m_err,mask,slp_m,slp_n)
      do j=js,je
         do i=is,ie
            if ( abs(phis(i,j)-gz0(i,j))/grav > 2. ) then
                 m_err(i,j) = mask(i,j)*(slp_m(i,j) - slp_n(i,j))*2.*grav/abs(phis(i,j)-gz0(i,j))
            else
                 m_err(i,j) = mask(i,j)*(slp_m(i,j) - slp_n(i,j))
            endif
         enddo
      enddo

      call rmse_bias(m_err, rms, bias, area)
      call corr(slp_m, slp_n, co, area)

      call prt_maxmin('SLP Error (mb)=', m_err, is, ie, js, je, 0, 1, 0.01)
      if(master) write(*,*) 'SLP (Pa): RMS=', rms, ' Bias=', bias
      if(master) write(*,*) 'SLP correlation=',co
  endif

  if ( nudge_winds ) then

       allocate (du_obs(is:ie,js:je,npz) )
       allocate (dv_obs(is:ie,js:je,npz) )

       du_obs = 0.
       dv_obs = 0.

! Compute tendencies:
     rdt = 1. / (tau_winds/factor + dt)
!$OMP parallel do default(none) shared(kstart,npz,kbot_winds,is,ie,js,je,du_obs,dv_obs, &
!$OMP                                  profile,u_obs,v_obs,ua,va,rdt,profile_Lower)
     do k=kstart, npz - kbot_winds
        do j=js,je
           do i=is,ie
              du_obs(i,j,k) = profile(k)*(u_obs(i,j,k)-ua(i,j,k))*rdt*profile_Lower(k)
              dv_obs(i,j,k) = profile(k)*(v_obs(i,j,k)-va(i,j,k))*rdt*profile_Lower(k)
           enddo
        enddo
     enddo
     if (nudge_option .ne. 0) call system_clock(c0)

     if ( nf_uv>0 .and. nudge_option.ne.2 .and. nudge_option.ne.4) call del2_uv(du_obs, dv_obs, del2_cd, npz, nf_uv, bd, npx, npy, gridstruct, domain)
     if (L_uv.ne.0 .and. nf_uv.gt.0 .and. (nudge_option.eq.2 .or. nudge_option.eq.4)) then
       if (do_nudge_now) then
         Lx_m = L_uv * 1000.0
         Ly_m = L_uv * 1000.0
         if (rf_debug_dump) call rf_debug_dump_slice('uv_before_rf',du_obs,npz, &
              agrid(is:ie,js:je,1),agrid(is:ie,js:je,2),dv_obs)
         call build_rf_alpha_from_L(gridstruct, Lx_m, Ly_m, bd, domain)
         call rf_uv(du_obs, dv_obs, kmd=npz, bd=bd, nf=nf_uv, npx=npx, npy=npy, &
              gridstruct=gridstruct, domain=domain)
         if (rf_debug_dump) call rf_debug_dump_slice('uv_after_rf',du_obs,npz, &
              agrid(is:ie,js:je,1),agrid(is:ie,js:je,2),dv_obs)
         du_stored=du_obs
         dv_stored=dv_obs
       else
         du_obs=du_stored
         dv_obs=dv_stored
       end if
     end if
     if (nudge_option .ne. 0) then
       call system_clock(c1)
       tcost = real(c1-c0)/real(rate)
       if (nudge_debug .and. is_master()) then
         write(*,'(a,f10.4,a)') 'RF2D cost = ', tcost, ' s'
       end if
     endif

! Apply U/V latitude taper AFTER spatial filtering.
if (nudge_option .ne. 0 .and. lat_nudge_taper_uv > 0.0) then
!$OMP parallel do default(none) shared(kstart,npz,kbot_winds,is,ie,js,je,du_obs,dv_obs,lat_weight_uv) private(i,j,k)
   do k=kstart, npz-kbot_winds
      do j=js,je
         do i=is,ie
            du_obs(i,j,k) = du_obs(i,j,k) * lat_weight_uv(i,j)
            dv_obs(i,j,k) = dv_obs(i,j,k) * lat_weight_uv(i,j)
         enddo
      enddo
   enddo
!$OMP end parallel do
endif

!$OMP parallel do default(none) shared(kstart,kbot_winds,npz,is,ie,js,je,du_obs,dv_obs, &
!$OMP                                  mask,ps_fac,u_dt,v_dt,ua,va,dt, climate_nudging)
     do k=kstart, npz - kbot_winds
       if ( climate_nudging ) then
         do j=js,je
           do i=is,ie
! Apply TC mask
             du_obs(i,j,k) = du_obs(i,j,k) * mask(i,j)
             dv_obs(i,j,k) = dv_obs(i,j,k) * mask(i,j)
             u_dt(i,j,k) = u_dt(i,j,k) + du_obs(i,j,k)
             v_dt(i,j,k) = v_dt(i,j,k) + dv_obs(i,j,k)
               ua(i,j,k) =   ua(i,j,k) + du_obs(i,j,k)*dt
               va(i,j,k) =   va(i,j,k) + dv_obs(i,j,k)*dt
           enddo
         enddo
       else
         if ( k==npz ) then
           do j=js,je
             do i=is,ie
               du_obs(i,j,k) = du_obs(i,j,k) * mask(i,j) * ps_fac(i,j)
               dv_obs(i,j,k) = dv_obs(i,j,k) * mask(i,j) * ps_fac(i,j)
!
               u_dt(i,j,k) = u_dt(i,j,k) + du_obs(i,j,k)
               v_dt(i,j,k) = v_dt(i,j,k) + dv_obs(i,j,k)
                 ua(i,j,k) =   ua(i,j,k) + du_obs(i,j,k)*dt
                 va(i,j,k) =   va(i,j,k) + dv_obs(i,j,k)*dt
             enddo
           enddo
         else
           do j=js,je
             do i=is,ie
! Apply TC mask
               du_obs(i,j,k) = du_obs(i,j,k) * mask(i,j)
               dv_obs(i,j,k) = dv_obs(i,j,k) * mask(i,j)
!
               u_dt(i,j,k) = u_dt(i,j,k) + du_obs(i,j,k)
               v_dt(i,j,k) = v_dt(i,j,k) + dv_obs(i,j,k)
                 ua(i,j,k) =   ua(i,j,k) + du_obs(i,j,k)*dt
                 va(i,j,k) =   va(i,j,k) + dv_obs(i,j,k)*dt
             enddo
           enddo
         endif
       endif ! climate_nudging
     enddo
  endif
!$OMP parallel do default(none) shared(is,ie,js,je,npz,t_dt)
  do k=1,npz
     do j=js,je
        do i=is,ie
           t_dt(i,j,k) = 0.
        enddo
     enddo
  enddo

  if ( nudge_virt ) then
       rdt = 1./(tau_virt/factor + dt)

!$OMP parallel do default(none) shared(is,ie,js,je,npz,kstart,kht,t_dt,prof_t,t_obs,zvir,profile_Lower, &
#ifdef MULTI_GASES
!$OMP                                  num_gas,                                           &
#endif
!$OMP                                  q,pt,rdt,ps_fac,climate_nudging)
     do k=kstart, kht
       if ( climate_nudging ) then
         do j=js,je
           do i=is,ie
             t_dt(i,j,k) = prof_t(k)*(t_obs(i,j,k)-pt(i,j,k))*rdt*profile_Lower(k)
           enddo
         enddo
       else
         if ( k==npz ) then
           do j=js,je
             do i=is,ie
#ifdef MULTI_GASES
                 t_dt(i,j,k) = prof_t(k)*(t_obs(i,j,k)/virq(q(i,j,k,1:num_gas))-pt(i,j,k))*rdt*ps_fac(i,j)*profile_Lower(k)
#else
                 t_dt(i,j,k) = prof_t(k)*(t_obs(i,j,k)/(1.+zvir*q(i,j,k,1))-pt(i,j,k))*rdt*ps_fac(i,j)*profile_Lower(k)
#endif
             enddo
           enddo
         else
           do j=js,je
             do i=is,ie
#ifdef MULTI_GASES
                 t_dt(i,j,k) = prof_t(k)*(t_obs(i,j,k)/virq(q(i,j,k,1:num_gas))-pt(i,j,k))*rdt*profile_Lower(k)
#else
                 t_dt(i,j,k) = prof_t(k)*(t_obs(i,j,k)/(1.+zvir*q(i,j,k,1))-pt(i,j,k))*rdt*profile_Lower(k)
#endif
             enddo
           enddo
         endif
       endif
     enddo
  endif

  if ( nudge_hght .and. kht<npz ) then     ! averaged (in log-p) temperature
       rdt = factor_nwp*1. / (tau_hght/factor + dt)
!$OMP parallel do default(none) shared(is,ie,js,je,npz,ak,h2,delp,kht,t_obs,q_obs,pt,zvir, &
#ifdef MULTI_GASES
!$OMP                                  num_gas,                &
#endif
!$OMP                                  q,rdt,ps_fac,mask,t_dt,climate_nudging) &
!$OMP                          private(pe2, peln )
       do j=js,je
           do i=is,ie
              pe2(i,1) = ak(1)
               h2(i,j) = 0.
           enddo
           do k=1, npz
              do i=is,ie
                 pe2(i,k+1) = pe2(i,k) + delp(i,j,k)
              enddo
           enddo
           do k=kht+1, npz+1
              do i=is,ie
                 peln(i,k) = log(pe2(i,k))
              enddo
           enddo

           do k=kht+1, npz
              do i=is,ie
! Difference between "Mean virtual tempearture" (pseudo height)
                if ( climate_nudging ) then
                  h2(i,j) = h2(i,j) + ((t_obs(i,j,k)*(1.+zvir*q_obs(i,j,k))-pt(i,j,k))*(1.+zvir*q(i,j,k,1))) &
                          *(peln(i,k+1)-peln(i,k))
                else
#ifdef MULTI_GASES
                  h2(i,j) = h2(i,j) + (t_obs(i,j,k)-pt(i,j,k)*virq(q(i,j,k,1:num_gas)))*(peln(i,k+1)-peln(i,k))
#else
                  h2(i,j) = h2(i,j) + (t_obs(i,j,k)-pt(i,j,k)*(1.+zvir*q(i,j,k,1)))*(peln(i,k+1)-peln(i,k))
#endif
                endif
              enddo
           enddo
           do i=is,ie
              h2(i,j) = h2(i,j) / (peln(i,npz+1)-peln(i,kht+1))
              h2(i,j) = rdt*ps_fac(i,j)*h2(i,j)*mask(i,j)
           enddo

           do k=kht+1, npz
              do i=is,ie
#ifdef MULTI_GASES
                 t_dt(i,j,k) = h2(i,j) / virq(q(i,j,k,1:num_gas))
#else
                 t_dt(i,j,k) = h2(i,j) / (1.+zvir*q(i,j,k,1))
#endif
              enddo
           enddo
       enddo   ! j-loop
       if(nudge_debug) call prt_maxmin('H2 increment=', h2, is, ie, js, je, 0, 1, dt)
  endif

  if ( nudge_virt .or. nudge_hght ) then
       if (nudge_option .ne. 0) call system_clock(c0)
       if ( nf_t>0 .and. nudge_option.ne.2 .and. nudge_option.ne.4) call del2_scalar(t_dt, del2_cd, npz, nf_t, bd, npx, npy, gridstruct, domain)
       if ( L_t.ne.0 .and. nf_t.gt.0 .and. (nudge_option.eq.2 .or. nudge_option.eq.4)) then
         if (do_nudge_now) then
           Lx_m = L_t * 1000.0
           Ly_m = L_t * 1000.0
           if (rf_debug_dump) call rf_debug_dump_slice('t_before_rf',t_dt,npz, &
                agrid(is:ie,js:je,1),agrid(is:ie,js:je,2))
           call build_rf_alpha_from_L(gridstruct, Lx_m, Ly_m, bd, domain)
           call rf2d_scalar(t_dt, kmd=npz, nf=nf_t, npx=npx, npy=npy, &
                   bd=bd, domain=domain, do_norm=.true.)
           if (rf_debug_dump) call rf_debug_dump_slice('t_after_rf',t_dt,npz, &
                agrid(is:ie,js:je,1),agrid(is:ie,js:je,2))
           dt_stored=t_dt
         else
           t_dt=dt_stored
         end if
       end if
! Apply temperature latitude taper AFTER spatial filtering.
! Skip the multiplication entirely when disabled.
     if (nudge_option .ne. 0 .and. lat_nudge_taper_t > 0.0) then
!$OMP parallel do default(none) shared(kstart,npz,is,ie,js,je,t_dt,lat_weight_t) private(i,j,k)
        do k=kstart, npz
           do j=js,je
              do i=is,ie
                 t_dt(i,j,k) = t_dt(i,j,k) * lat_weight_t(i,j)
              enddo
           enddo
        enddo
!$OMP end parallel do
     endif

     if (nudge_option .ne. 0) then
       call system_clock(c1)
       tcost = real(c1-c0)/real(rate)
       if (nudge_debug .and. is_master()) then
         write(*,'(a,f10.4,a)') 'RF2D cost2 = ', tcost, ' s'
       end if
     endif

!$OMP parallel do default(none) shared(kstart,npz,is,ie,js,je,pt,t_dt,dt,mask)
       do k=kstart, npz
          do j=js,je
             do i=is,ie
                pt(i,j,k) = pt(i,j,k) + t_dt(i,j,k)*dt*mask(i,j)
            enddo
         enddo
       enddo
  endif

  !!q_dt(:,:,:) = 0.
  if ( nudge_q ) then
       rdt = factor_nwp*1./(tau_q/factor + dt)
!$OMP parallel do default(none) shared(kstart,npz,kbot_q,is,ie,js,je,press,p_wvp,nwat, &
!$OMP                                  q,delp,q_dt,prof_q,q_min,q_obs,rdt,mask,dt,climate_nudging)
     do k=kstart, npz - kbot_q
        if ( press(k) > p_wvp ) then
            do iq=2,nwat
               do j=js,je
                  do i=is,ie
                     q(i,j,k,iq) = q(i,j,k,iq)*delp(i,j,k)
                  enddo
               enddo
            enddo
! Specific humidity:
            do j=js,je
               do i=is,ie
                  delp(i,j,k) = delp(i,j,k)*(1.-q(i,j,k,1))
                  if ( climate_nudging ) then
                    q_dt(i,j,k,1) = q_dt(i,j,k,1)+prof_q(k)*(max(q_min,q_obs(i,j,k))-q(i,j,k,1))*rdt
                    q(i,j,k,1) = q(i,j,k,1) + prof_q(k)*(max(q_min,q_obs(i,j,k))-q(i,j,k,1))*rdt*dt
                  else
                    q_dt(i,j,k,1) = prof_q(k)*(max(q_min,q_obs(i,j,k))-q(i,j,k,1))*rdt*mask(i,j)
                    q(i,j,k,1) = q(i,j,k,1) + q_dt(i,j,k,1)*dt
                  endif
                  delp(i,j,k) = delp(i,j,k)/(1.-q(i,j,k,1))
               enddo
            enddo
            do iq=2,nwat
               do j=js,je
                  do i=is,ie
                     q(i,j,k,iq) = q(i,j,k,iq)/delp(i,j,k)
                  enddo
               enddo
            enddo
        endif
     enddo
  endif

  ps_dt(is:ie,js:je) = 0.0

  deallocate ( t_obs )
  deallocate ( q_obs )
  deallocate ( ps_obs )

  if ( breed_srf_w .and. nudge_winds )   &
  call breed_srf_winds(Time, dt, npz, u_obs, v_obs, ak, bk, ps, phis, delp, ua, va, u_dt, v_dt, pt, q, nwat, zvir, gridstruct)

  if ( nudge_debug) then
       call prt_maxmin('T increment=', t_dt,   is, ie, js, je, 0, npz, dt)
       if ( nudge_winds ) then
            call prt_maxmin('U increment=', du_obs, is, ie, js, je, 0, npz, dt)
            call prt_maxmin('V increment=', dv_obs, is, ie, js, je, 0, npz, dt)
       endif
  endif

  if (nudge_option .ne. 0 .and. nf_skip_timestep > 1) skip_counter=skip_counter-1
  if ( nudge_winds .or. nudge_option .ne. 0 ) then
     deallocate ( u_obs )
     deallocate ( v_obs )
  endif
  if ( nudge_winds ) then
     deallocate (du_obs)
     deallocate (dv_obs)
  endif

  nullify(agrid)
  nullify(area)
  nullify(rarea)

  nullify(vlon)
  nullify(vlat)
  nullify(sina_u)
  nullify(sina_v)
  nullify(sin_sg)

  nullify(dx)
  nullify(dy)
  nullify(rdxc)
  nullify(rdyc)

  nullify(da_min)

  nullify(sw_corner)
  nullify(se_corner)
  nullify(nw_corner)
  nullify(ne_corner)

 end  subroutine fv_nwp_nudge


 subroutine ps_nudging(dt, factor, factor_nwp, npz, ak, bk, ps_obs, mask, tm, ps, phis, delp, ua, va, pt, nwat, q, bd, npx, npy, gridstruct, domain)
! Input
      real, intent(in):: dt, factor, factor_nwp
      integer, intent(in):: npz, nwat, npx, npy
      real, intent(in), dimension(npz+1):: ak, bk
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(in):: phis(isd:ied,jsd:jed)
      real, intent(in), dimension(is:ie,js:je):: ps_obs, mask, tm
      type(fv_grid_type), intent(IN), target :: gridstruct
      type(domain2d), intent(INOUT) :: domain
! Input/Output
      real, intent(inout), dimension(isd:ied,jsd:jed):: ps
      real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp, pt, ua, va
      real, intent(inout):: q(isd:ied,jsd:jed,npz,nwat)
! local
      real, dimension(is:ie,js:je):: ps_dt
      integer, parameter:: kmax = 100
      real:: pn0(kmax+1), pk0(kmax+1), pe0(kmax+1)
      real, dimension(is:ie,npz+1):: pe2, peln
      real:: pst, dbk, pt0, rdt, bias
      integer i, j, k, iq

      real, pointer, dimension(:,:) :: area
#ifdef MULTI_GASES
      real :: kappax(km)
      real :: pkx
#endif

      area => gridstruct%area

! Adjust interpolated ps to model terrain
    if ( kmax < km ) call mpp_error(FATAL,'==> KMAX must be larger than km')

    do j=js,je
       do i=is,ie
#ifdef MULTI_GASES
       do k=1,km
          kappax(k)= virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))
       enddo
#endif
       do k=1, km+1
          pe0(k) =  ak0(k) + bk0(k)*ps_obs(i,j)
          pk0(k) = (ak0(k) + bk0(k)*ps_obs(i,j))**kappa
       enddo
      if( phis(i,j)>gz0(i,j) ) then
          do k=km,1,-1
             if( phis(i,j)<gz3(i,j,k) .and. phis(i,j) >= gz3(i,j,k+1) ) then
                 pst = pk0(k) + (pk0(k+1)-pk0(k))*(gz3(i,j,k)-phis(i,j))/(gz3(i,j,k)-gz3(i,j,k+1))
                 go to 666
             endif
          enddo
      else
          pn0(km  ) = log(ak0(km) + bk0(km)*ps_obs(i,j))
          pn0(km+1) = log(ps_obs(i,j))
! Extrapolation into the ground using only the lowest layer potential temp
#ifdef MULTI_GASES
           pkx = (pk0(km+1)-pk0(km))/(kappa*(pn0(km+1)-pn0(km)))
           pt0 = tm(i,j)/exp(kappax(km)*log(pkx))
           pkx = exp((kappax(km)-1.0)*log(pkx))
           pst = pk0(km+1) + (gz0(i,j)-phis(i,j))/(cp_air*pt0*pkx)
#else
           pt0 = tm(i,j)/(pk0(km+1)-pk0(km))*(kappa*(pn0(km+1)-pn0(km)))
           pst = pk0(km+1) + (gz0(i,j)-phis(i,j))/(cp_air*pt0)
#endif
      endif
#ifdef MULTI_GASES
      k=km
666   ps_dt(i,j) = pst**(1./(kappa*kappax(k))) - ps(i,j)
#else
666   ps_dt(i,j) = pst**(1./kappa) - ps(i,j)
#endif
      enddo   ! i-loop
      enddo   ! j-loop

      if( nf_ps>0 ) call del2_scalar(ps_dt, del2_cd, 1, nf_ps, bd, npx, npy, gridstruct, domain)

      do j=js,je
         do i=is,ie
! Cap the errors:
            ps_dt(i,j) = sign ( min(10.E2, abs(ps_dt(i,j))), ps_dt(i,j) )
            ps_dt(i,j) = mask(i,j)*ps_dt(i,j) * max( 0., 1.-abs(gz0(i,j)-phis(i,j))/(grav*500.) )
         enddo
      enddo

      if( do_ps_bias ) call ps_bias_correction( ps_dt, is, ie, js, je, isd, ied, jsd, jed, area)

#ifdef CONSV_HGHT
! Convert tracer moist mixing ratio to mass
       do iq=1,nwat
       do k=1,npz
          do j=js,je
             do i=is,ie
                q(i,j,k,iq) = q(i,j,k,iq) * delp(i,j,k)
             enddo
          enddo
       enddo
       enddo

       do k=1,npz
          do j=js,je
             do i=is,ie
                ua(i,j,k) = ua(i,j,k) * delp(i,j,k)
                va(i,j,k) = va(i,j,k) * delp(i,j,k)
             enddo
          enddo
       enddo

      do j=js,je
         do i=is,ie
            pe2(i,1) = ak(1)
            peln(i,1) = log(pe2(i,1))
         enddo
         do k=2,npz+1
             do i=is,ie
                 pe2(i,k) = pe2(i,k-1) + delp(i,j,k-1)
                peln(i,k) = log(pe2(i,k))
             enddo
         enddo
         do k=1,npz
             do i=is,ie
                pt(i,j,k) = pt(i,j,k)*(peln(i,k+1)-peln(i,k))
             enddo
         enddo
      enddo
#endif

! Update ps:
      do j=js,je
         do i=is,ie
            ps(i,j) = ak(1)
         enddo
      enddo

      rdt = factor_nwp*dt / (tau_ps/factor + dt)
      do k=1,npz
         dbk = rdt*(bk(k+1) - bk(k))
         do j=js,je
            do i=is,ie
               delp(i,j,k) = delp(i,j,k) + dbk*ps_dt(i,j)
                   ps(i,j) = delp(i,j,k) + ps(i,j)
            enddo
        enddo
     enddo

#ifdef CONSV_HGHT
      do j=js,je
         do k=2,npz+1
             do i=is,ie
                pe2(i,k) = pe2(i,k-1) + delp(i,j,k-1)
                peln(i,k) = log(pe2(i,k))
             enddo
         enddo
         do k=1,npz
             do i=is,ie
                pt(i,j,k) = pt(i,j,k)/(peln(i,k+1)-peln(i,k))
             enddo
          enddo
      enddo

! Convert tracer density back to moist mixing ratios
       do iq=1,nwat
       do k=1,npz
          do j=js,je
             do i=is,ie
                q(i,j,k,iq) = q(i,j,k,iq) / delp(i,j,k)
             enddo
          enddo
       enddo
       enddo

       do k=1,npz
          do j=js,je
             do i=is,ie
                ua(i,j,k) = ua(i,j,k) / delp(i,j,k)
                va(i,j,k) = va(i,j,k) / delp(i,j,k)
             enddo
          enddo
       enddo
#endif

 end subroutine ps_nudging

 subroutine ps_bias_correction ( ps_dt, is, ie, js, je, isd, ied, jsd, jed, area )
 integer, intent(IN) :: is,  ie,  js, je
 integer, intent(IN) :: isd, ied, jsd,jed
 real, intent(inout):: ps_dt(is:ie,js:je)
 real, intent(IN), dimension(isd:ied,jsd:jed) :: area
!
 real:: esl, total_area
 real:: bias, psum
 integer:: i, j

 total_area = 4.*pi*radius**2
 esl = 0.01       ! Pascal

 bias = g0_sum(ps_dt, is, ie, js, je, 1, .true., isd, ied, jsd, jed, area)

 if ( abs(bias) < esl ) then
      if(master .and. nudge_debug) write(*,*) 'Very small PS bias=', -bias, ' No bais adjustment'
      return
 else
      if(master .and. nudge_debug) write(*,*) 'Significant PS bias=', -bias
 endif

 if ( bias > 0. ) then
     psum = 0.
     do j=js,je
        do i=is,ie
           if ( ps_dt(i,j) > 0. ) then
                psum = psum + area(i,j)
           endif
        enddo
     enddo

     call mp_reduce_sum( psum )
     bias = bias * total_area / psum

     do j=js,je
        do i=is,ie
           if ( ps_dt(i,j) > 0.0 ) then
                ps_dt(i,j) = max(0.0, ps_dt(i,j) - bias)
           endif
        enddo
     enddo
 else
     psum = 0.
     do j=js,je
        do i=is,ie
           if ( ps_dt(i,j) < 0. ) then
                psum = psum + area(i,j)
           endif
        enddo
     enddo

     call mp_reduce_sum( psum )
     bias = bias * total_area / psum

     do j=js,je
        do i=is,ie
           if ( ps_dt(i,j) < 0.0 ) then
                ps_dt(i,j) = min(0.0, ps_dt(i,j) - bias)
           endif
        enddo
     enddo
 endif

 end subroutine ps_bias_correction

!>@brief Fast version of global sum that is reproduced
!! if the result is rounded to 4-byte
 real function g0_sum(p, ifirst, ilast, jfirst, jlast, mode, reproduce, isd, ied, jsd, jed, area)
! Fast version of global sum; reproduced if result rounded to 4-byte
      integer, intent(IN) :: ifirst, ilast
      integer, intent(IN) :: jfirst, jlast
      integer, intent(IN) :: mode  !< if ==1 divided by global area
      logical, intent(IN) :: reproduce
      real,    intent(IN) :: p(ifirst:ilast,jfirst:jlast)      !< field to be summed
      integer, intent(IN) :: isd, ied, jsd, jed
      real,    intent(IN) :: area(isd:ied,jsd:jed)

      integer :: i,j
      real gsum

      gsum = 0.
      do j=jfirst,jlast
         do i=ifirst,ilast
            gsum = gsum + p(i,j)*area(i,j)
         enddo
      enddo
!     call mpp_sum(gsum)    ! does this work?
      call mp_reduce_sum(gsum)

      if ( mode==1 ) then
           g0_sum = gsum / (4.*pi*radius**2)
      else
           g0_sum = gsum
      endif

      if ( reproduce ) g0_sum = real(g0_sum, 4) ! convert to 4-byte real

 end function g0_sum


 subroutine compute_slp(isc, iec, jsc, jec, tm, ps, gz, slp)
 integer, intent(in):: isc, iec, jsc, jec
 real, intent(in), dimension(isc:iec,jsc:jec):: tm, ps, gz
! tm: virtual temperature required as input
 real, intent(out):: slp(isc:iec,jsc:jec)
 integer:: i,j

    do j=jsc,jec
       do i=isc,iec
          slp(i,j) = ps(i,j) * exp( gz(i,j)/(rdgas*(tm(i,j) + 3.25E-3*gz(i,j)/grav)) )
       enddo
    enddo

 end subroutine compute_slp


 subroutine get_obs(Time, dt, zvir, ak, bk, ps, ts, ps_obs, delp, pt, nwat, q, u_obs, v_obs, t_obs, q_obs,  &
                    phis, ua, va, u_dt, v_dt, npx, npy, npz, factor, factor_nwp, mask, ptop, bd, gridstruct, domain, tile_num)
  type(time_type), intent(in):: Time
  integer,         intent(in):: npz, nwat, npx, npy, tile_num
  real,            intent(in):: zvir, ptop
  real,            intent(in):: dt, factor, factor_nwp
  real, intent(in), dimension(npz+1):: ak, bk
  type(fv_grid_bounds_type), intent(IN) :: bd
  real, intent(in), dimension(isd:ied,jsd:jed):: phis
  real, intent(in), dimension(is:ie,js:je):: mask
  real, intent(inout), dimension(isd:ied,jsd:jed):: ps
  real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp, pt, ua, va, u_dt, v_dt
  real, intent(inout)::q(isd:ied,jsd:jed,npz,*)
  real, intent(out), dimension(is:ie,js:je):: ts, ps_obs
  real, intent(out), dimension(is:ie,js:je,npz):: u_obs, v_obs, t_obs, q_obs
  type(fv_grid_type), intent(IN) :: gridstruct
  type(domain2d), intent(INOUT) :: domain
! local:
  real::  tm(is:ie,js:je)
  real(KIND=4), allocatable,dimension(:,:,:):: ut, vt, wt
  real, allocatable,dimension(:,:,:):: uu, vv
  integer :: seconds, days
  integer :: i,j,k
  real :: alpha, beta
  character(len=256) :: fname_tile
  character(len=8)   :: c_tile

  call get_time (time, seconds, days)

  ! New input paths do not provide analyzed surface temperature.
  ! Preserve the legacy nudge_option=0 behavior; define a fallback only for
  ! the new nudging options.
  if ( nudge_option .ne. 0 ) ts(:,:) = pt(is:ie,js:je,npz)

  if ( do_adiabatic_init  ) goto 333
  seconds = seconds - nint(dt)

! Data must be "time_interval" (hr) apart; keep two time levels in memory

  no_obs = .false.
  analysis_time = .false.

  if ( mod(seconds, time_interval) == 0 ) then

    if ( nfile == nfile_total ) then
         no_obs = .true.
#ifndef DYCORE_SOLO
         forecast_mode = .true.
#endif
         if(print_end_nudge)  then
            print_end_nudge = .false.
            if (master) write(*,*) '*** L-S nudging Ended at', days, seconds
         endif
         return              ! free-running mode
    endif

      ps_dat(:,:,1) = ps_dat(:,:,2)
      if ( nudge_winds ) then
         u_dat(:,:,:,1) = u_dat(:,:,:,2)
         v_dat(:,:,:,1) = v_dat(:,:,:,2)
      endif
      t_dat(:,:,:,1) = t_dat(:,:,:,2)
      q_dat(:,:,:,1) = q_dat(:,:,:,2)
      if ( nudge_option.ne.0 ) then
         dp_dat(:,:,:,1) = dp_dat(:,:,:,2)
      endif

!---------------
! Read next data
!---------------
      nfile = nfile + 1
      if ( nudge_option.eq.0 ) then
        call get_ncep_analysis ( ps_dat(:,:,2), u_dat(:,:,:,2), v_dat(:,:,:,2),    &
                              t_dat(:,:,:,2), q_dat(:,:,:,2), zvir,  &
                              ts, nfile, file_names(nfile), bd )
      else if ( grib2_input ) then
        call get_grib2_analysis( ps_dat(:,:,2), u_dat(:,:,:,2), v_dat(:,:,:,2), &
                              t_dat(:,:,:,2), q_dat(:,:,:,2), dp_dat(:,:,:,2), &
                              nfile, file_names(nfile), domain, zvir )
      else if ( nudge_option.eq.3 .or. nudge_option.eq.4) then
        write(c_tile, '(I0)') tile_num
        fname_tile = trim(file_names(nfile)) // ".tile" // trim(adjustl(c_tile))
        call get_regional_analysis( ps_dat(:,:,2), u_dat(:,:,:,2), v_dat(:,:,:,2), &
                            t_dat(:,:,:,2), q_dat(:,:,:,2), dp_dat(:,:,:,2), &
                            nfile, fname_tile, bd, gridstruct, ak, bk, &
                            npx, npy, -1 , domain, zvir)
      else if ( nudge_option.eq.1 .or. nudge_option.eq.2 ) then
        call get_regional_analysis( ps_dat(:,:,2), u_dat(:,:,:,2), v_dat(:,:,:,2), &
                            t_dat(:,:,:,2), q_dat(:,:,:,2), dp_dat(:,:,:,2), &
                            nfile, file_names(nfile), bd, gridstruct, ak, bk, &
                            npx, npy, -1 , domain, zvir)
      else
        if(master)  write(*,*) 'Nudging regional option not available, WIP'
      end if
      time_nudge = dt
  else
      time_nudge = time_nudge + dt
  endif

!--------------------
! Time interpolation:
!--------------------

  beta = time_nudge / real(time_interval)

  if ( beta < 0. .or. beta >  (1.+1.E-7) ) then
       if (master) write(*,*) 'Nudging: beta=', beta
       call mpp_error(FATAL,'==> Error from get_obs:data out of range')
  endif

333  continue

  if ( do_adiabatic_init ) then
       beta = 1.; alpha = 0.
      if( nudge_debug .and. master) write(*,*) 'Doing final adiabatic initialization/nudging'
  else
      alpha = 1. - beta
      if( abs(alpha)<1.e-7 )  analysis_time = .true.
  endif

! Warning: ps_data are not adjusted for the differences in terrain yet.
! Pressure-level GRIB2 input has no analysis surface pressure.  Do not let the
! internal placeholder ps_dat contaminate ps_fac or diagnostics; use the actual
! current model surface pressure on this path.  nudge_ps is prohibited for GRIB2.
  if ( grib2_input ) then
     ps_obs(:,:) = ps(is:ie,js:je)
  else
     ps_obs(:,:) = alpha*ps_dat(:,:,1) + beta*ps_dat(:,:,2)
  endif

!---------------------------------
!*** nudge & update ps & delp here
!---------------------------------
  if (nudge_ps) then

      allocate ( wt(is:ie,js:je,km) )
      wt(:,:,:) = alpha*t_dat(:,:,:,1) + beta*t_dat(:,:,:,2)
! Needs gz3 for ps_nudging
      call get_int_hght(npz, ak, bk, ps(is:ie,js:je), delp, ps_obs(is:ie,js:je), wt)
      do j=js,je
         do i=is,ie
            tm(i,j)  = wt(i,j,km)
        enddo
      enddo
      deallocate ( wt )

      allocate ( uu(isd:ied,jsd:jed,npz) )
      allocate ( vv(isd:ied,jsd:jed,npz) )
      uu = ua
      vv = va
      call ps_nudging(dt, factor, factor_nwp, npz, ak, bk, ps_obs, mask, tm, ps, phis, delp, uu, vv, pt, nwat, q, bd, npx, npy, gridstruct, domain)
      do k=1,npz
      do j=js,je
         do i=is,ie
            u_dt(i,j,k) = u_dt(i,j,k) + (uu(i,j,k) - ua(i,j,k)) / dt
            v_dt(i,j,k) = v_dt(i,j,k) + (vv(i,j,k) - va(i,j,k)) / dt
        enddo
      enddo
      enddo
      deallocate (uu )
      deallocate (vv )
  endif

  allocate ( ut(is:ie,js:je,npz) )
  allocate ( vt(is:ie,js:je,npz) )

  if ( nudge_winds ) then

     if ( nudge_option.eq.0 ) then
       call remap_uv(npz, ak,  bk, ps(is:ie,js:je), delp,  ut,     vt,   &
                     km, ps_dat(is:ie,js:je,1),  u_dat(:,:,:,1), v_dat(:,:,:,1), ptop )
     else if ( grib2_input ) then
       call remap_uv_grib2(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, grib2_plev_pa, u_dat(:,:,:,1), v_dat(:,:,:,1), ptop)
     else
       call remap_uv_regional(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, dp_dat(:,:,:,1), u_dat(:,:,:,1), v_dat(:,:,:,1), ptop )
     endif

       u_obs(:,:,:) = alpha*ut(:,:,:)
       v_obs(:,:,:) = alpha*vt(:,:,:)

     if ( nudge_option.eq.0 ) then
       call remap_uv(npz, ak, bk, ps(is:ie,js:je), delp,   ut,      vt,   &
                     km, ps_dat(is:ie,js:je,2),  u_dat(:,:,:,2), v_dat(:,:,:,2), ptop )
     else if ( grib2_input ) then
       call remap_uv_grib2(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, grib2_plev_pa, u_dat(:,:,:,2), v_dat(:,:,:,2), ptop)
     else
       call remap_uv_regional(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, dp_dat(:,:,:,2), u_dat(:,:,:,2), v_dat(:,:,:,2), ptop )
     endif

       u_obs(:,:,:) = u_obs(:,:,:) + beta*ut(:,:,:)
       v_obs(:,:,:) = v_obs(:,:,:) + beta*vt(:,:,:)
  endif
     if ( nudge_option.eq.0 ) then
       call remap_tq(npz, ak, bk, ps(is:ie,js:je), delp,  ut,  vt,  &
                     km,  ps_dat(is:ie,js:je,1),  t_dat(:,:,:,1), q_dat(:,:,:,1), zvir, ptop)
     else if ( grib2_input ) then
       call remap_tq_grib2(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, grib2_plev_pa, t_dat(:,:,:,1), q_dat(:,:,:,1), zvir, ptop)
     else
       call remap_tq_regional(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, dp_dat(:,:,:,1), t_dat(:,:,:,1), q_dat(:,:,:,1), zvir, ptop)
     endif

       t_obs(:,:,:) = alpha*ut(:,:,:)
       q_obs(:,:,:) = alpha*vt(:,:,:)

     if ( nudge_option.eq.0 ) then
       call remap_tq(npz, ak, bk, ps(is:ie,js:je), delp,  ut,  vt,  &
                     km,  ps_dat(is:ie,js:je,2),  t_dat(:,:,:,2), q_dat(:,:,:,2), zvir, ptop)
     else if ( grib2_input ) then
       call remap_tq_grib2(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, grib2_plev_pa, t_dat(:,:,:,2), q_dat(:,:,:,2), zvir, ptop)
     else
       call remap_tq_regional(npz, ak, bk, ps(is:ie,js:je), delp, ut, vt, &
                     km, dp_dat(:,:,:,2), t_dat(:,:,:,2), q_dat(:,:,:,2), zvir, ptop)
     endif

       t_obs(:,:,:) = t_obs(:,:,:) + beta*ut(:,:,:)
       q_obs(:,:,:) = q_obs(:,:,:) + beta*vt(:,:,:)

  deallocate ( ut )
  deallocate ( vt )

 end subroutine get_obs

!>@brief The subroutine 'fv_nwp_nudge_init' initializes the nudging module.
!>@details This module opens analysis files and computes remapping coefficients.
 subroutine fv_nwp_nudge_init(time, axes, npz, zvir, ak, bk, ts, phis, gridstruct, ks, npx, npy, neststruct, bd, domain, tile_num)
 character(len=17) :: mod_name = 'fv_nudge'
  type(time_type), intent(in):: time
  integer,         intent(in):: axes(4)
  integer,  intent(in):: npz           !< vertical dimension
  integer,  intent(in):: tile_num        !< Tile Number
  real,     intent(in):: zvir
  type(fv_grid_bounds_type), intent(IN) :: bd
  real, intent(in), dimension(bd%isd:bd%ied,bd%jsd:bd%jed):: phis
  real, intent(in), dimension(npz+1):: ak, bk
  real, intent(out), dimension(bd%is:bd%ie,bd%js:bd%je):: ts
  type(fv_grid_type), target :: gridstruct
  integer,  intent(in) :: ks, npx, npy
  type(fv_nest_type) :: neststruct
  type(domain2d), intent(INOUT), target :: domain

  real :: missing_value = -1.e10
  logical found
  integer tsize(4)
  integer :: i, j, j1, f_unit, unit, io, ierr, nt, k
  integer :: ncid
  character(len=256) :: fname_tile
  character(len=8)   :: c_tile

! ---> h1g, read NCEP analysis data list, 2012-10-22
  integer            :: input_fname_unit
  character(len=128) :: fname_tmp
! <--- h1g, 2012-10-22

  real, pointer, dimension(:,:,:) :: agrid

  is  = bd%is
  ie  = bd%ie
  js  = bd%js
  je  = bd%je

  isd = bd%isd
  ied = bd%ied
  jsd = bd%jsd
  jed = bd%jed


   agrid => gridstruct%agrid

   master = is_master()
   do_adiabatic_init = .false.
   deg2rad = pi/180.
   rad2deg = 180./pi

   if (neststruct%nested) then
!!! Assumes no grid stretching
      grid_size = 1.E7/(neststruct%npx_global*neststruct%refinement_of_global)
   else
      grid_size = 1.E7/real(npx-1)         ! mean grid size
   endif

   do nt=1,nfile_max
      file_names(nt) = "No_File_specified"
   enddo

   track_file_name = "No_File_specified"

   ! Default GRIB2 pressure-level table.  The namelist may override the
   ! active level count and the leading grib2_plev_pa entries.
   grib2_nplev = grib2_nplev_default
   grib2_plev_pa = 0.0
   grib2_plev_pa(1:grib2_nplev_default) = grib2_default_plev_pa

    if( file_exists( 'input.nml' ) ) then
        read( input_nml_file, nml = fv_nwp_nudge_nml, iostat = io )
        ierr = check_nml_error(io,'fv_nwp_nudge_nml')
    end if
    if (nudge_option .ne. 0) then
       no_obs = .false.
       analysis_time = .false.
       time_nudge = 0.0
       skip_counter = 0
       nudge_tile_num = tile_num
    endif
    if (nudge_option .ne. 0) then
       if (nf_skip_grid < 0) then
          if (master) write(*,*) 'WARNING: nf_skip_grid < 0; using 0'
          nf_skip_grid = 0
       else if (nf_skip_grid > 3) then
          if (master) write(*,*) 'WARNING: nf_skip_grid > 3; be careful when nudge_option = 2 with small subdomains.'
       endif
    endif

    if (grib2_input) then
       if (nudge_option.eq.0) then
          call mpp_error(FATAL,'grib2_input is intended for nudge_option > 0')
       endif
       if (grib2_nplev < 2 .or. grib2_nplev > grib2_nplev_max) then
          call mpp_error(FATAL,'grib2_nplev must be between 2 and grib2_nplev_max')
       endif
       do k=1,grib2_nplev
          if (grib2_plev_pa(k) <= 0.0) then
             call mpp_error(FATAL,'grib2_plev_pa active levels must be positive')
          endif
       enddo
       do k=2,grib2_nplev
          if (grib2_plev_pa(k) <= grib2_plev_pa(k-1)) then
             call mpp_error(FATAL,'grib2_plev_pa must increase downward')
          endif
       enddo
    endif
    call write_version_number ( 'FV_NUDGE_MOD', version )
    if ( master ) then
         f_unit=stdlog()
         write( f_unit, nml = fv_nwp_nudge_nml )
         write(*,*) 'NWP nudging initialized.'
    endif
! ---> h1g, specify the NCEP analysis data for nudging, 2012-10-22
    if ( trim(input_fname_list) .ne. "" ) then
      input_fname_unit = get_unit()
      open( input_fname_unit, file = input_fname_list )
       nt = 0
       io = 0

       do while ( io .eq. 0 )
        read ( input_fname_unit, '(a)', iostat = io, end = 101 ) fname_tmp
        if( trim(fname_tmp) .ne. "" ) then    ! escape any empty record
          if ( trim(fname_tmp) == trim(analysis_file_last) ) then
            nt = nt + 1
            file_names(nt) = trim(fname_tmp)
            if(master .and. nudge_debug) write(*,*) 'From NCEP file list, last file: ', nt, file_names(nt)
            nt = 0
            goto 101  ! read last analysis data and then close file
          endif ! trim(fname_tmp) == trim(analysis_file_last)

          if ( trim(analysis_file_first) == "" ) then
            nt = nt + 1
            file_names(nt) = trim(fname_tmp)
            if(master .and. nudge_debug) then
              if( nt .eq. 1 ) then
               write(*,*) 'From NCEP file list, first file: ', nt, file_names(nt),trim(analysis_file_first)
              else
               write(*,*) 'From NCEP file list: ', nt, file_names(nt)
              endif ! nt .eq. 1
            endif ! master .and. nudge_debug
          else
            if ( trim(fname_tmp) == trim(analysis_file_first) .or. nt > 0 ) then
              nt = nt + 1
              file_names(nt) = trim(fname_tmp)
              if(master .and. nudge_debug) then
                if( nt .eq. 1 ) then
                  write(*,*) 'From NCEP file list, first file: ', nt,  file_names(nt),trim(analysis_file_first)
                else
                  write(*,*) 'From NCEP file list: ', nt, file_names(nt)
                endif !  nt .eq. 1
              endif  ! master .and. nudge_debug
            endif  ! trim(fname_tmp) == trim(analysis_file_first) .or. nt > 0
          endif  ! trim(analysis_file_first) == ""
        endif  ! trim(fname_tmp) .ne. ""
       end do !  io .eq. 0
101   close( input_fname_unit )
    endif
! <--- h1g, 2012-10-22

    do nt=1,nfile_max
      if ( file_names(nt) == "No_File_specified" ) then
           nfile_total = nt - 1
           if(master) write(*,*) 'Total of NCEP files specified=', nfile_total
           exit
      endif
    enddo

    id_ht_err = register_diag_field ( mod_name, 'ht_error', axes(1:2), time,        &
                'height_error', 'DEG K', missing_value=missing_value )


! Initialize remapping coefficients:

!   call field_size(file_names(1), 'T', tsize, field_found=found)
!   if ( found ) then
!        im = tsize(1); jm = tsize(2); km = tsize(3)
!        if(master)  write(*,*) 'NCEP analysis dimensions:', tsize
!   else
!        call mpp_error(FATAL,'==> Error from get_ncep_analysis: T field not found')
!   endif

    if ( nudge_option.eq.0 ) then
      call open_ncfile( file_names(1), ncid )        ! open the file
      call get_ncdim1( ncid, 'lon', im )
      call get_ncdim1( ncid, 'lat', jm )
      call get_ncdim1( ncid, 'lev', km )
      if(master)  write(*,*) 'NCEP analysis dimensions:', im, jm, km

      allocate ( s2c(is:ie,js:je,4) )
      allocate ( id1(is:ie,js:je) )
      allocate ( id2(is:ie,js:je) )
      allocate ( jdc(is:ie,js:je) )

      allocate (  lon(im) )
      allocate (  lat(jm) )

      call _GET_VAR1 (ncid, 'lon', im, lon )
      call _GET_VAR1 (ncid, 'lat', jm, lat )

!   Convert to radian
      do i=1,im
         lon(i) = lon(i) * deg2rad ! lon(1) = 0.
      enddo
      do j=1,jm
         lat(j) = lat(j) * deg2rad
      enddo

      allocate ( ak0(km+1) )
      allocate ( bk0(km+1) )

      call _GET_VAR1 (ncid, 'hyai', km+1, ak0, found )
      if ( .not. found )  ak0(:) = 0.

      call _GET_VAR1 (ncid, 'hybi', km+1, bk0 )
      call close_ncfile( ncid )

!   Note: definition of NCEP hybrid is p(k) = a(k)*1.E5 + b(k)*ps
      if ( .not. using_merra2) then
        ! This is not needed for MERRA2 data
        ak0(:) = ak0(:) * 1.E5
      endif
!   Limiter to prevent NAN at top during remapping
      if ( bk0(1) < 1.E-9 ) ak0(1) = max(1.e-9, ak0(1))

     if ( master ) then
        do k=1,npz
           write(*,*) k, 0.5*(ak(k)+ak(k+1))+0.5*(bk(k)+bk(k+1))*1.E5,  'del-B=', bk(k+1)-bk(k)
        enddo
     endif

     if ( k_breed==0 ) k_breed = max(1, ks)

     call slp_obs_init

!  -----------------------------------------------------------
!   Initialize lat-lon to Cubed bi-linear interpolation coeff:
!  -----------------------------------------------------------
      call remap_coef(agrid)

!   Find bounding latitudes:
      jbeg = jm-1;         jend = 2
      do j=js,je
         do i=is,ie
              j1 = jdc(i,j)
            jbeg = min(jbeg, j1)
            jend = max(jend, j1+1)
         enddo
      enddo

      allocate ( gz0(is:ie,js:je) )
      allocate ( gz3(is:ie,js:je,km+1) )
      allocate (ps_dat(is:ie,js:je,2) )
      allocate ( u_dat(is:ie,js:je,km,2) )
      allocate ( v_dat(is:ie,js:je,km,2) )
      allocate ( t_dat(is:ie,js:je,km,2) )
      allocate ( q_dat(is:ie,js:je,km,2) )


!   Get first dataset
      nt = 2
      nfile = 1
      call get_ncep_analysis ( ps_dat(:,:,nt), u_dat(:,:,:,nt), v_dat(:,:,:,nt),     &
                              t_dat(:,:,:,nt), q_dat(:,:,:,nt), zvir,   &
                              ts, nfile, file_names(nfile), bd )
    else
      if ( grib2_input ) then
        if ( nudge_ps ) then
          call mpp_error(FATAL,'Pressure-level GRIB2 nudging does not provide surface pressure; set nudge_ps=.false.')
        endif
        ! Determine the regular lat-lon source grid from the GRIB2 metadata
        ! and build the existing bilinear lat-lon -> FV3 A-grid weights.
        km = grib2_nplev
        call grib2_init_grid(file_names(1), agrid, domain)
      else
        if (nudge_option.eq.3 .or. nudge_option.eq.4) then
          write(c_tile, '(I0)') tile_num
          fname_tile = trim(file_names(1)) // ".tile" // trim(adjustl(c_tile))
          call open_ncfile( fname_tile, ncid )
        else
          call open_ncfile( file_names(1), ncid )
        end if
        call get_ncdim1( ncid, 'lon', im )
        call get_ncdim1( ncid, 'lat', jm )
        call get_ncdim1( ncid, 'lev', km )
        call close_ncfile( ncid )
      endif

      allocate ( gz0(is:ie,js:je) )
      allocate (ps_dat(is:ie,js:je,2) )
      ps_dat=0.0
      allocate (dp_dat(is:ie,js:je,km,2) )
      allocate ( u_dat(is:ie,js:je,km,2) )
      allocate ( v_dat(is:ie,js:je,km,2) )
      allocate ( t_dat(is:ie,js:je,km,2) )
      allocate ( q_dat(is:ie,js:je,km,2) )
      if ( nudge_ps ) allocate ( gz3(is:ie,js:je,km+1) )

      nt = 2
      nfile = 1
      if ( grib2_input ) then
        call get_grib2_analysis( ps_dat(:,:,nt), u_dat(:,:,:,nt), v_dat(:,:,:,nt), &
                            t_dat(:,:,:,nt), q_dat(:,:,:,nt), dp_dat(:,:,:,nt), &
                            nfile, file_names(nfile), domain, zvir )
      else if (nudge_option.eq.3 .or. nudge_option.eq.4) then
        write(c_tile, '(I0)') tile_num
        fname_tile = trim(file_names(nfile)) // ".tile" // trim(adjustl(c_tile))
        call get_regional_analysis( ps_dat(:,:,nt), u_dat(:,:,:,nt), v_dat(:,:,:,nt), &
                            t_dat(:,:,:,nt), q_dat(:,:,:,nt), dp_dat(:,:,:,nt), &
                            nfile, fname_tile, bd, gridstruct, ak, bk, &
                            npx, npy, -1 , domain, zvir )
      else
        call get_regional_analysis( ps_dat(:,:,nt), u_dat(:,:,:,nt), v_dat(:,:,:,nt), &
                            t_dat(:,:,:,nt), q_dat(:,:,:,nt), dp_dat(:,:,:,nt), &
                            nfile, file_names(nfile), bd, gridstruct, ak, bk, &
                            npx, npy, -1 , domain, zvir )
      end if

      ! Retained for the existing regional get_int_hght/remap interfaces.
      ! The pressure-level GRIB2 path itself does not use ak0/bk0.
      allocate ( ak0(km+1) )
      allocate ( bk0(km+1) )
      ak0 = 0.0
      bk0 = 0.0
      ak0(1:min(km+1,npz+1)) = ak(1:min(km+1,npz+1))
      bk0(1:min(km+1,npz+1)) = bk(1:min(km+1,npz+1))
    end if !End nudge options

    module_is_initialized = .true.
    nullify(agrid)


 end subroutine fv_nwp_nudge_init


 subroutine get_ncep_analysis ( ps, u, v, t, q, zvir, ts, nfile, fname, bd )
  real,     intent(in):: zvir
  character(len=128), intent(in):: fname
  integer,  intent(inout):: nfile
!
  type(fv_grid_bounds_type), intent(IN) :: bd
  real, intent(out), dimension(is:ie,js:je):: ts, ps
  real(KIND=4), intent(out), dimension(is:ie,js:je,km):: u, v, t, q
! local:
  real(kind=4), allocatable:: wk0(:,:), wk1(:,:), wk2(:,:), wk3(:,:,:)
  real tmean
  integer:: i, j, k, npt
  integer:: i1, i2, j1, ncid
  logical found
  logical:: read_ts = .true.
  logical:: land_ts = .false.

  integer:: status, var3id    ! h1g, 2016-08-10
#include <netcdf.inc>


  if( .not. file_exists(fname) ) then
     call mpp_error(FATAL,'==> Error from get_ncep_analysis: file not found: '//fname)
  else
     call open_ncfile( fname, ncid )        ! open the file
     if(master) write(*,*) 'Reading NCEP anlysis file:', fname
  endif

  if ( climate_nudging ) read_ts =.false.
  if ( read_ts ) then       ! read skin temperature; could be used for SST
       allocate ( wk1(im,jm) )

       call get_var3_r4( ncid, 'TS', 1,im, 1,jm, 1,1, wk1 )
!      if ( master ) write(*,*) 'Done reading NCEP TS data'

      if ( .not. land_ts ) then
           allocate ( wk0(im,jm) )
! Read NCEP ORO (1; land; 0: ocean; 2: sea_ice)

! ---> h1g, read either 'ORO' or 'LAND', 2016-08-10
        status = nf_inq_varid (ncid,  'ORO', var3id)
        if (status .eq. NF_NOERR)  then
          call get_var3_r4( ncid, 'ORO', 1,im, 1,jm, 1,1, wk0 )

        else !there is no 'ORO'
          status = nf_inq_varid (ncid,  'LAND', var3id)
          if (status .eq. NF_NOERR)  then
            call get_var3_r4( ncid, 'LAND', 1,im, 1,jm, 1,1, wk0 )
          else
            call mpp_error(FATAL,'Neither ORO nor LAND exists in re-analysis data')
          endif
        endif
! <--- h1g, 2016-08-10

           do j=1,jm
              tmean = 0.
              npt = 0
              do i=1,im
                 if( abs(wk0(i,j)-1.) > 0.99 ) then ! non-land points
                     tmean = tmean + wk1(i,j)
                     npt = npt + 1
                 endif
              enddo
!-------------------------------------------------------
! Replace TS over interior land with zonal mean SST/Ice
!-------------------------------------------------------
              if ( npt /= 0 ) then
                   tmean= tmean / real(npt)
                   do i=1,im
                      if( abs(wk0(i,j)-1.) <= 0.99 ) then   ! land points
                          if ( i==1 ) then
                               i1 = im;     i2 = 2
                          elseif ( i==im ) then
                               i1 = im-1;   i2 = 1
                          else
                               i1 = i-1;    i2 = i+1
                          endif
                          if ( abs(wk0(i2,j)-1.)>0.99 ) then     ! east side has priority
                               wk1(i,j) = wk1(i2,j)
                          elseif ( abs(wk0(i1,j)-1.)>0.99 ) then ! west side
                               wk1(i,j) = wk1(i1,j)
                          else
                               wk1(i,j) = tmean
                          endif
                      endif
                   enddo
              endif
           enddo
           deallocate ( wk0 )
      endif   ! land_ts

      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            ts(i,j) = s2c(i,j,1)*wk1(i1,j1  ) + s2c(i,j,2)*wk1(i2,j1  ) +  &
                      s2c(i,j,3)*wk1(i2,j1+1) + s2c(i,j,4)*wk1(i1,j1+1)
         enddo
      enddo
      call prt_maxmin('SST_model', ts, is, ie, js, je, 0, 1, 1.)

#ifndef DYCORE_SOLO
! Perform interp to FMS SST format/grid
      call ncep2fms( wk1 )
      if(master) call pmaxmin( 'SST_ncep', real(sst_ncep), i_sst, j_sst, 1.)
!     if(nfile/=1 .and. master) call pmaxmin( 'SST_anom', sst_anom, i_sst, j_sst, 1.)
#endif
       deallocate ( wk1 )
       if (master) write(*,*) 'Done processing NCEP SST'

  endif     ! read_ts

!----------------------------------
! remap surface pressure and height:
!----------------------------------

     allocate ( wk2(im,jbeg:jend) )
     call get_var3_r4( ncid, 'PS', 1,im, jbeg,jend, 1,1, wk2 )

     do j=js,je
        do i=is,ie
           i1 = id1(i,j)
           i2 = id2(i,j)
           j1 = jdc(i,j)
           ps(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                     s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
        enddo
     enddo


! ---> h1g, read either 'PHIS' or 'PHI', 2016-08-10
     status = nf_inq_varid (ncid,  'PHIS', var3id)
     if (status .eq. NF_NOERR)  then
       call get_var3_r4( ncid, 'PHIS', 1,im, jbeg,jend, 1,1, wk2 )

     else !there is no 'PHIS'
       status = nf_inq_varid (ncid,  'PHI', var3id)
       if (status .eq. NF_NOERR)  then
         call get_var3_r4( ncid, 'PHI', 1,im, jbeg,jend, 1,1, wk2 )
         wk2 = wk2 * grav  ! convert unit from geopotential meter (m) to geopotential height (m2/s2)
       else
         call mpp_error(FATAL,'Neither PHIS nor PHI exists in re-analysis data')
       endif
     endif
! <--- h1g, 2016-08-10


     do j=js,je
        do i=is,ie
           i1 = id1(i,j)
           i2 = id2(i,j)
           j1 = jdc(i,j)
           gz0(i,j) = s2c(i,j,1)*wk2(i1,j1  ) + s2c(i,j,2)*wk2(i2,j1  ) +  &
                      s2c(i,j,3)*wk2(i2,j1+1) + s2c(i,j,4)*wk2(i1,j1+1)
        enddo
     enddo
     call prt_maxmin('ZS_ncep', gz0, is, ie, js, je, 0, 1, 1./grav)
     deallocate ( wk2 )


   allocate ( wk3(1:im, jbeg:jend, 1:km) )

! Winds:
   if ( nudge_winds ) then

      call get_var3_r4( ncid, 'U', 1,im, jbeg,jend, 1,km, wk3 )

      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            u(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

      call get_var3_r4( ncid, 'V', 1,im, jbeg,jend, 1,km, wk3 )

      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            v(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

   endif

! Read in tracers: only sphum at this point
      call get_var3_r4( ncid, 'Q', 1,im, jbeg,jend, 1,km , wk3 )

      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            q(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

      call get_var3_r4( ncid, 'T', 1,im, jbeg,jend, 1,km , wk3 )
      call close_ncfile ( ncid )

      do k=1,km
      do j=js,je
         do i=is,ie
            i1 = id1(i,j)
            i2 = id2(i,j)
            j1 = jdc(i,j)
            t(i,j,k) = s2c(i,j,1)*wk3(i1,j1  ,k) + s2c(i,j,2)*wk3(i2,j1  ,k) +  &
                       s2c(i,j,3)*wk3(i2,j1+1,k) + s2c(i,j,4)*wk3(i1,j1+1,k)
         enddo
      enddo
      enddo

      if ( .not. climate_nudging) then
      if ( .not. T_is_Tv ) then
      do k=1,km
      do j=js,je
         do i=is,ie
! The field T in Larry H.'s post processing of NCEP analysis is actually virtual temperature
! before Dec 1, 2005
! Convert t to virtual temperature:
            t(i,j,k) = t(i,j,k)*(1.+zvir*q(i,j,k))
         enddo
      enddo
      enddo
      endif
      endif

!  endif

   deallocate ( wk3 )

! nfile = nfile + 1

 end subroutine get_ncep_analysis

subroutine get_regional_analysis ( ps, u, v, t, q, dp, nfile, fname, &
                                   bd, gridstruct, ak, bk,           &
                                   npx, npy, grid_type, domain, zvir)
!-----------------------------------------------------------------------
!  Reads a regional FV3 analysis on the D-grid (u_s,v_s,u_w,v_w),
!  rotates all winds into Earth-relative (east/north) using full
!  spherical geometry, converts to A-grid, and passes vertical levels
!  "as-is" to the nudging system.
!
!  No vertical interpolation is performed here.
!-----------------------------------------------------------------------
  use fv_arrays_mod,     only: fv_grid_type, fv_grid_bounds_type
  use fv_grid_utils_mod, only: mid_pt_sphere, get_unit_vect2, get_latlon_vector, inner_prod, cubed_to_latlon
  use mpp_domains_mod,   only: domain2d
  use fms2_io_mod,       only: file_exists
  use mpp_mod,           only: mpp_error, FATAL, mpp_pe, mpp_root_pe
  use sim_nc_mod,        only: open_ncfile, close_ncfile, get_ncdim1, get_var3_r4
!  use fv_nesting_mod,    only: d2a_setup
  implicit none

  !---------------- arguments ----------------
  character(len=128),          intent(in)    :: fname
  integer,                     intent(inout) :: nfile
  type(fv_grid_bounds_type),   intent(in)    :: bd
  type(fv_grid_type),          intent(in), target :: gridstruct
  real,                        intent(in)    :: ak(:), bk(:)
  integer,                     intent(in)    :: npx,npy,grid_type
  type(domain2d), intent(inout) :: domain
  real,                        intent(in)    :: zvir

  real,          intent(out) :: ps(is:ie, js:je)
  real(kind=4),  intent(out) :: u (is:ie, js:je, km)
  real(kind=4),  intent(out) :: v (is:ie, js:je, km)
  real(kind=4),  intent(out) :: t (is:ie, js:je, km)
  real(kind=4),  intent(out) :: q (is:ie, js:je, km)
  real(kind=4),  intent(out) :: dp(is:ie, js:je, km)

  !---------------- locals ----------------
  integer :: i, j, k, lev, levp, ncid, kmax
  real, pointer :: grid(:,:,:), cosa_s(:,:), rsin2(:,:)
  real(kind=4) :: glon(is:ie,js:je),glat(is:ie,js:je)
  ! Input analysis fields read by get_var3_r4 must remain explicit real(4).
  real(kind=4), allocatable :: ps_in(:,:)
  real(kind=4), allocatable :: t_in(:,:,:), q_in(:,:,:), delp_in(:,:,:)
  real(kind=4), allocatable :: u_s(:,:,:), v_s(:,:,:)   ! S/N edges
  real(kind=4), allocatable :: u_w(:,:,:), v_w(:,:,:)   ! W/E edges

  ! Rotated D-grid -> Earth-relative
  real, allocatable :: ud(:,:,:), vd(:,:,:)

  ! A-grid Earth-relative
  real, allocatable :: u_a(:,:,:), v_a(:,:,:)
  real, allocatable :: u_d(:,:,:), v_d(:,:,:)
  real, allocatable :: u_e(:,:,:), v_e(:,:,:)

  ! geometry vectors
  real(kind=R_GRID) :: p1(2), p2(2), p3(2)
  real(kind=R_GRID) :: e1(3), e2(3)
  real(kind=R_GRID) :: ex(3), ey(3), ex_c(3), ey_c(3)
  real(kind=R_GRID) :: V3(3)

  !--------- extract FV3 grid ------------
  grid => gridstruct%grid   ! (isd:ied+1, jsd:jed+1, 2) lon/lat in radians
  cosa_s => gridstruct%cosa_s
  rsin2  => gridstruct%rsin2
  !========= Read file ==========
  if (.not. file_exists(fname)) then
     call mpp_error(FATAL,'get_regional_analysis: file not found '//trim(fname))
  endif

  call open_ncfile(fname,ncid)

  call get_ncdim1(ncid,'lev', lev)
  call get_ncdim1(ncid,'levp',levp)

  !========= Allocate input ===========
  allocate(ps_in   (is:ie, js:je))
  allocate(t_in    (is:ie, js:je, lev))
  allocate(q_in    (is:ie, js:je, lev))
  allocate(delp_in (is:ie, js:je, lev))

  allocate(u_s(is:ie,   js:je+1, lev))
  allocate(v_s(is:ie,   js:je+1, lev))
  allocate(u_w(is:ie+1, js:je,   lev))
  allocate(v_w(is:ie+1, js:je,   lev))

  ! Rotated edges
  allocate(ud(is:ie+1, js:je+1, lev))
  allocate(vd(is:ie+1, js:je+1, lev))

  ! A-grid
  allocate(u_e(is:ie, js:je, lev))
  allocate(v_e(is:ie, js:je, lev))

  !========= Read variables ==========
  call get_var3_r4(ncid,'geolon',     is,ie,js,je,1,1,glon)
  call get_var3_r4(ncid,'geolat',     is,ie,js,je,1,1,glat)
  call get_var3_r4(ncid,'ps',     is,ie,js,je,1,1,ps_in)
  ps(is:ie,js:je) = ps_in
  call get_var3_r4(ncid,'t',      is,ie,js,je,1,lev, t_in)
  call get_var3_r4(ncid,'sphum',  is,ie,js,je,1,lev, q_in)
  call get_var3_r4(ncid,'delp',   is,ie,js,je,1,lev, delp_in)

  call get_var3_r4(ncid,'u_s', is,ie,   js,je+1,1,lev,  u_s)
  call get_var3_r4(ncid,'v_s', is,ie,   js,je+1,1,lev,  v_s)
  call get_var3_r4(ncid,'u_w', is,ie+1, js,je,  1,lev,  u_w)
  call get_var3_r4(ncid,'v_w', is,ie+1, js,je,  1,lev,  v_w)

  call close_ncfile(ncid)

!====================================================================
! 1. Rotate S/N and W/E edge winds into Earth-relative east/north
!====================================================================
  do k = 1, lev
     ! South/North edges: ud(i,j,k) along local east direction
     do j = js, je+1
        do i = is, ie
           p1(:) = grid(i  , j, 1:2)
           p2(:) = grid(i+1, j, 1:2)
           call mid_pt_sphere(p1, p2, p3)
           call get_unit_vect2(p1, p2, e1)        ! tangent along edge
           call get_latlon_vector(p3, ex, ey)     ! earth east/north
           ud(i,j,k) = real( &
                u_s(i,j,k)*inner_prod(e1,ex) + &
                v_s(i,j,k)*inner_prod(e1,ey), &
                kind=R_GRID )
        enddo
     enddo
     ! West/East edges: vd(i,j,k) along local north direction
     do j = js, je
        do i = is, ie+1
           p1(:) = grid(i, j  , 1:2)
           p2(:) = grid(i, j+1, 1:2)
           call mid_pt_sphere(p1, p2, p3)
           call get_unit_vect2(p1, p2, e2)        ! tangent along edge
           call get_latlon_vector(p3, ex, ey)
           vd(i,j,k) = real( &
                u_w(i,j,k)*inner_prod(e2,ex) + &
                v_w(i,j,k)*inner_prod(e2,ey), &
                kind=R_GRID )
        enddo
     enddo
  enddo
!====================================================================
! 2. A-grid wind: average grid-relative D-grid points
!====================================================================
  allocate( u_d(isd:ied,   jsd:jed+1, lev) )
  allocate( v_d(isd:ied+1, jsd:jed,   lev) )
  allocate( u_a(isd:ied, jsd:jed, lev) )
  allocate( v_a(isd:ied, jsd:jed, lev) )
  u_d = 0.0
  v_d = 0.0
  ! Embed our local subdomain ud,vd into u_d,v_d
  do k = 1, lev
     do j = js, je+1
        do i = is, ie
           u_d(i,j,k) = ud(i,j,k)
        enddo
     enddo
     do j = js, je
        do i = is, ie+1
           v_d(i,j,k) = vd(i,j,k)
        enddo
     enddo
  enddo
  call cubed_to_latlon( u_d, v_d, u_a, v_a, &
                        gridstruct, npx, npy, lev,  &
                        0, grid_type, domain, gridstruct%bounded_domain, 4 , bd ) !NOTE FOR CHECK LATER: Is c2l_ord always 4 in HAFS regional?
  do k = 1, lev
     do j = js, je
        do i = is, ie
           u_e(i,j,k) = u_a(i,j,k)
           v_e(i,j,k) = v_a(i,j,k)
        enddo
     enddo
  end do
!====================================================================
! 3. Vertical copy to output (no remap: nudging handles that)
!====================================================================
  kmax = min(km, lev)

  u  = 0.0
  v  = 0.0
  t  = 0.0
  q  = 0.0
  dp = 0.0

  do k = 1, kmax
     do j = js, je
        do i = is, ie
           u(i,j,k)  = u_e(i,j,k)
           v(i,j,k)  = v_e(i,j,k)
           q(i,j,k)  = q_in(i,j,k)
           if ( climate_nudging ) then
             if ( T_is_Tv ) then
               t(i,j,k) = t_in(i,j,k) / (1. + zvir*q_in(i,j,k))
             else
               t(i,j,k) = t_in(i,j,k)
             endif
           else
             if ( .not. T_is_Tv ) then
               t(i,j,k) = t_in(i,j,k) * (1. + zvir*q_in(i,j,k))
             else
               t(i,j,k) = t_in(i,j,k)
             endif
           endif
           dp(i,j,k) = delp_in(i,j,k)
        enddo
     enddo
  enddo

  if (km > kmax) then
     do k = kmax+1, km
        do j = js, je
           do i = is, ie
              u(i,j,k)  = u_e(i,j,kmax)
              v(i,j,k)  = v_e(i,j,kmax)
              q(i,j,k)  = q_in(i,j,kmax)
              if ( climate_nudging ) then
                if ( T_is_Tv ) then
                  t(i,j,k) = t_in(i,j,kmax) / (1. + zvir*q_in(i,j,kmax))
                else
                  t(i,j,k) = t_in(i,j,kmax)
                endif
              else
                if ( .not. T_is_Tv ) then
                  t(i,j,k) = t_in(i,j,kmax) * (1. + zvir*q_in(i,j,kmax))
                else
                  t(i,j,k) = t_in(i,j,kmax)
                endif
              endif
              dp(i,j,k) = delp_in(i,j,kmax)
           enddo
        enddo
     enddo
  endif

end subroutine get_regional_analysis


!=======================================================================
! Direct pressure-level GRIB2 support for nudge_option > 0.
!
! The GRIB2 files are expected to contain regular lat-lon, instantaneous
! isobaric fields using GRIB2 PDT 4.0:
!   TMP  : discipline 0, category 0, parameter 0
!   SPFH : discipline 0, category 1, parameter 0
!   UGRD : discipline 0, category 2, parameter 2
!   VGRD : discipline 0, category 2, parameter 3
!   level type 100 (isobaric surface), level value in Pa.
!
! SPFH is read even when nudge_q=.false. because the existing nudging
! code normally carries temperature as virtual temperature.  Thus TMP
! (true temperature) is converted to Tv before storage except when
! climate_nudging expects true temperature.
!=======================================================================

subroutine grib2_init_grid(fname, agrid, domain)
  use mpp_mod,         only: mpp_pe, mpp_send, mpp_recv, mpp_sync_self, mpp_error, FATAL
  use mpp_domains_mod, only: domain2d
  implicit none
  character(len=*), intent(in) :: fname
  real, intent(in) :: agrid(isd:ied,jsd:jed,2)
  type(domain2d), intent(inout) :: domain

  integer :: mype, root, npelist, n, ierr_g2, lugb
  integer, allocatable :: pelist(:)
  real, allocatable :: sample(:,:), sample_valid(:,:), lon_deg(:), lat_deg(:)
  real :: dlon, lon_tol, lat_tol

  mype = mpp_pe()
  call get_rf_domain_pelist(domain,pelist,npelist,root)

  if (mype == root) then
    lugb = 91
    call baopenr(lugb,trim(fname),ierr_g2)
    if (ierr_g2 /= 0) then
      call mpp_error(FATAL,'grib2_init_grid: cannot open '//trim(fname))
    endif

    call grib2_read_field_open(lugb,0,0,nint(grib2_plev_pa(1)), &
                                     sample,sample_valid,lon_deg,lat_deg,.true.,.true.,.true.)

    im = size(sample,1)
    jm = size(sample,2)

    call baclose(lugb,ierr_g2)
    call gf_finalize(ierr_g2)
  endif

  if (mype == root) then
    do n=1,npelist
      if (pelist(n) /= root) then
        call mpp_send(im,plen=1,to_pe=pelist(n))
        call mpp_send(jm,plen=1,to_pe=pelist(n))
      endif
    enddo
  else
    call mpp_recv(im,glen=1,from_pe=root,block=.true.)
    call mpp_recv(jm,glen=1,from_pe=root,block=.true.)
  endif
  call mpp_sync_self()

  if (allocated(lon)) deallocate(lon)
  if (allocated(lat)) deallocate(lat)
  allocate(lon(im),lat(jm))

  if (mype == root) then
    lon(:) = lon_deg(:)*deg2rad
    lat(:) = lat_deg(:)*deg2rad
  endif
  call rf_pelist_broadcast_real(lon(1),im,root,pelist,npelist)
  call rf_pelist_broadcast_real(lat(1),jm,root,pelist,npelist)

  ! remap_coef uses a cyclic longitude closure.  Do not silently apply that
  ! assumption to a regional GRIB2 grid.  The current direct-GRIB2 path is
  ! therefore restricted to global regular longitude coverage, which is the
  ! AIGFS use case.  Latitude must also cover every local FV3 target point.
  if (im < 2 .or. jm < 2) then
    call mpp_error(FATAL,'GRIB2 source grid must have at least 2x2 points')
  endif
  dlon=lon(2)-lon(1)
  lon_tol=max(1.e-8,1.e-3*abs(dlon))
  if (abs((lon(im)-lon(1)+dlon)-2.*pi) > lon_tol) then
    call mpp_error(FATAL,'GRIB2 direct nudging currently requires global longitude coverage')
  endif
  lat_tol=1.e-5  ! radians: r4-safe roundoff tolerance; still strict coverage
  if (minval(agrid(is:ie,js:je,2)) < lat(1)-lat_tol .or. &
      maxval(agrid(is:ie,js:je,2)) > lat(jm)+lat_tol) then
    call mpp_error(FATAL,'GRIB2 latitude coverage does not contain the local FV3 domain')
  endif

  if (allocated(s2c)) deallocate(s2c)
  if (allocated(id1)) deallocate(id1)
  if (allocated(id2)) deallocate(id2)
  if (allocated(jdc)) deallocate(jdc)
  allocate(s2c(is:ie,js:je,4))
  allocate(id1(is:ie,js:je))
  allocate(id2(is:ie,js:je))
  allocate(jdc(is:ie,js:je))

  call remap_coef(agrid)

  if (mype == root .and. nudge_debug) then
    write(*,'(a,2i8)') 'GRIB2 source grid im,jm = ',im,jm
    write(*,'(a,2f12.4)') 'GRIB2 longitude range = ', &
                          minval(lon_deg),maxval(lon_deg)
    write(*,'(a,2f12.4)') 'GRIB2 latitude range  = ', &
                          minval(lat_deg),maxval(lat_deg)
  endif

  if (mype == root) then
    if (allocated(sample)) deallocate(sample)
    if (allocated(sample_valid)) deallocate(sample_valid)
    if (allocated(lon_deg)) deallocate(lon_deg)
    if (allocated(lat_deg)) deallocate(lat_deg)
  endif
  deallocate(pelist)
end subroutine grib2_init_grid

subroutine get_grib2_analysis(ps,u,v,t,q,dp,nfile,fname,domain,zvir)
  use mpp_mod,         only: mpp_pe, mpp_error, FATAL
  use mpp_domains_mod, only: domain2d
  implicit none
  character(len=*), intent(in) :: fname
  integer, intent(inout) :: nfile
  type(domain2d), intent(inout) :: domain
  real, intent(in) :: zvir

  real, intent(out) :: ps(is:ie,js:je)
  real(kind=4), intent(out) :: u(is:ie,js:je,km)
  real(kind=4), intent(out) :: v(is:ie,js:je,km)
  real(kind=4), intent(out) :: t(is:ie,js:je,km)
  real(kind=4), intent(out) :: q(is:ie,js:je,km)
  real(kind=4), intent(out) :: dp(is:ie,js:je,km)

  integer :: i,j,k,i1,i2,j1
  integer :: mype,root,npelist,lugb,ierr_g2
  integer, allocatable :: pelist(:)
  real, allocatable :: src(:,:), src_valid(:,:), lon_dummy(:), lat_dummy(:)
  logical, allocatable :: u_ok(:,:,:), v_ok(:,:,:), t_ok(:,:,:), q_ok(:,:,:)
  real :: c1,c2,c3,c4,wsum

  if (.not. file_exists(fname)) then
    call mpp_error(FATAL,'get_grib2_analysis: file not found '//trim(fname))
  endif
  if (km /= grib2_nplev) then
    call mpp_error(FATAL,'get_grib2_analysis: km must equal grib2_nplev')
  endif

  mype = mpp_pe()
  call get_rf_domain_pelist(domain,pelist,npelist,root)

  allocate(src(im,jm),src_valid(im,jm))
  allocate(u_ok(is:ie,js:je,km),v_ok(is:ie,js:je,km))
  allocate(t_ok(is:ie,js:je,km),q_ok(is:ie,js:je,km))
  src = 0.0
  src_valid = 1.0
  u_ok=.false.; v_ok=.false.; t_ok=.false.; q_ok=.false.

  if (mype == root) then
    lugb = 91
    call baopenr(lugb,trim(fname),ierr_g2)
    if (ierr_g2 /= 0) then
      call mpp_error(FATAL,'get_grib2_analysis: cannot open '//trim(fname))
    endif
    if (nudge_debug) write(*,*) 'Reading pressure-level GRIB2: ',trim(fname)
  endif

  ! Pressure-level GRIB2 does not provide source surface pressure in the
  ! required input set.  ps is only a storage placeholder here; get_obs
  ! explicitly uses the actual model ps and sets ps_fac=1 for grib2_input.
  ps = 100000.0
  dp = 0.0
  u  = 0.0
  v  = 0.0
  t  = 0.0
  q  = 0.0

  do k=1,km
    ! UGRD: category 2, parameter 2.
    if (mype == root) then
      call grib2_read_field_open(lugb,2,2,nint(grib2_plev_pa(k)), &
                                       src,src_valid,lon_dummy,lat_dummy,.false.,k==1,.true.)
    endif
    call rf_pelist_broadcast_real(src(1,1),im*jm,root,pelist,npelist)
    call rf_pelist_broadcast_real(src_valid(1,1),im*jm,root,pelist,npelist)
    do j=js,je
      do i=is,ie
        i1=id1(i,j); i2=id2(i,j); j1=jdc(i,j)
        c1=s2c(i,j,1); c2=s2c(i,j,2); c3=s2c(i,j,3); c4=s2c(i,j,4)
        wsum=c1*src_valid(i1,j1)+c2*src_valid(i2,j1)+ &
             c3*src_valid(i2,j1+1)+c4*src_valid(i1,j1+1)
        if (wsum > 1.e-12) then
          u(i,j,k)=(c1*src_valid(i1,j1)*src(i1,j1)+ &
                    c2*src_valid(i2,j1)*src(i2,j1)+ &
                    c3*src_valid(i2,j1+1)*src(i2,j1+1)+ &
                    c4*src_valid(i1,j1+1)*src(i1,j1+1))/wsum
          u_ok(i,j,k)=.true.
        endif
      enddo
    enddo

    ! VGRD: category 2, parameter 3.
    if (mype == root) then
      call grib2_read_field_open(lugb,2,3,nint(grib2_plev_pa(k)), &
                                       src,src_valid,lon_dummy,lat_dummy,.false.,.false.,.true.)
    endif
    call rf_pelist_broadcast_real(src(1,1),im*jm,root,pelist,npelist)
    call rf_pelist_broadcast_real(src_valid(1,1),im*jm,root,pelist,npelist)
    do j=js,je
      do i=is,ie
        i1=id1(i,j); i2=id2(i,j); j1=jdc(i,j)
        c1=s2c(i,j,1); c2=s2c(i,j,2); c3=s2c(i,j,3); c4=s2c(i,j,4)
        wsum=c1*src_valid(i1,j1)+c2*src_valid(i2,j1)+ &
             c3*src_valid(i2,j1+1)+c4*src_valid(i1,j1+1)
        if (wsum > 1.e-12) then
          v(i,j,k)=(c1*src_valid(i1,j1)*src(i1,j1)+ &
                    c2*src_valid(i2,j1)*src(i2,j1)+ &
                    c3*src_valid(i2,j1+1)*src(i2,j1+1)+ &
                    c4*src_valid(i1,j1+1)*src(i1,j1+1))/wsum
          v_ok(i,j,k)=.true.
        endif
      enddo
    enddo

    ! TMP: category 0, parameter 0.  GRIB2 TMP is true temperature.
    if (mype == root) then
      call grib2_read_field_open(lugb,0,0,nint(grib2_plev_pa(k)), &
                                       src,src_valid,lon_dummy,lat_dummy,.false.,.false.,.true.)
    endif
    call rf_pelist_broadcast_real(src(1,1),im*jm,root,pelist,npelist)
    call rf_pelist_broadcast_real(src_valid(1,1),im*jm,root,pelist,npelist)
    do j=js,je
      do i=is,ie
        i1=id1(i,j); i2=id2(i,j); j1=jdc(i,j)
        c1=s2c(i,j,1); c2=s2c(i,j,2); c3=s2c(i,j,3); c4=s2c(i,j,4)
        wsum=c1*src_valid(i1,j1)+c2*src_valid(i2,j1)+ &
             c3*src_valid(i2,j1+1)+c4*src_valid(i1,j1+1)
        if (wsum > 1.e-12) then
          t(i,j,k)=(c1*src_valid(i1,j1)*src(i1,j1)+ &
                    c2*src_valid(i2,j1)*src(i2,j1)+ &
                    c3*src_valid(i2,j1+1)*src(i2,j1+1)+ &
                    c4*src_valid(i1,j1+1)*src(i1,j1+1))/wsum
          t_ok(i,j,k)=.true.
        endif
      enddo
    enddo

    ! SPFH: category 1, parameter 0.  It is needed for Tv in the normal
    ! nudging path and for q_obs in climate height nudging, even when q itself
    ! is not nudged.  Individual pressure levels are optional and are filled
    ! vertically below; an entirely missing q profile remains a fatal error.
    if (mype == root) then
      call grib2_read_field_open(lugb,1,0,nint(grib2_plev_pa(k)), &
                                       src,src_valid,lon_dummy,lat_dummy,.false.,.false.,.false.)
    endif
    call rf_pelist_broadcast_real(src(1,1),im*jm,root,pelist,npelist)
    call rf_pelist_broadcast_real(src_valid(1,1),im*jm,root,pelist,npelist)
    do j=js,je
      do i=is,ie
        i1=id1(i,j); i2=id2(i,j); j1=jdc(i,j)
        c1=s2c(i,j,1); c2=s2c(i,j,2); c3=s2c(i,j,3); c4=s2c(i,j,4)
        wsum=c1*src_valid(i1,j1)+c2*src_valid(i2,j1)+ &
             c3*src_valid(i2,j1+1)+c4*src_valid(i1,j1+1)
        if (wsum > 1.e-12) then
          q(i,j,k)=(c1*src_valid(i1,j1)*src(i1,j1)+ &
                    c2*src_valid(i2,j1)*src(i2,j1)+ &
                    c3*src_valid(i2,j1+1)*src(i2,j1+1)+ &
                    c4*src_valid(i1,j1+1)*src(i1,j1+1))/wsum
          q_ok(i,j,k)=.true.
        endif
      enddo
    enddo
  enddo

  ! NCEPLIBS-g2 expands bitmap-missing values to zero.  Do not let those
  ! zeros enter the nudging fields.  Horizontal interpolation above is
  ! renormalized over valid corners, and any pressure level still missing at
  ! a target point is filled vertically from the nearest valid profile data.
  call fill_grib2_vertical_missing(u,u_ok,km,grib2_plev_pa,'UGRD')
  call fill_grib2_vertical_missing(v,v_ok,km,grib2_plev_pa,'VGRD')
  call fill_grib2_vertical_missing(t,t_ok,km,grib2_plev_pa,'TMP')
  call fill_grib2_vertical_missing(q,q_ok,km,grib2_plev_pa,'SPFH')
  q=max(0.0,q)

  ! The normal nudging path carries t_obs as virtual temperature.  The
  ! climate_nudging branch expects true temperature, so leave TMP unchanged
  ! there and retain q separately for the height-nudging calculation.
  if ( .not. climate_nudging ) then
    do k=1,km
      do j=js,je
        do i=is,ie
          t(i,j,k)=t(i,j,k)*(1.0+zvir*q(i,j,k))
        enddo
      enddo
    enddo
  endif

  if (mype == root) then
    call baclose(lugb,ierr_g2)
    call gf_finalize(ierr_g2)
    if (nudge_debug) then
      write(*,'(a,2es16.6)') 'GRIB2 U min/max = ',minval(u),maxval(u)
      write(*,'(a,2es16.6)') 'GRIB2 V min/max = ',minval(v),maxval(v)
      if (climate_nudging) then
        write(*,'(a,2es16.6)') 'GRIB2 T min/max = ',minval(t),maxval(t)
      else
        write(*,'(a,2es16.6)') 'GRIB2 Tv min/max = ',minval(t),maxval(t)
      endif
    endif
  endif

  deallocate(src,src_valid,pelist)
  deallocate(u_ok,v_ok,t_ok,q_ok)
end subroutine get_grib2_analysis

subroutine fill_grib2_vertical_missing(field,valid,kmd,plev_pa,vname)
  use mpp_mod, only: mpp_error, FATAL, mpp_pe
  implicit none
  integer, intent(in) :: kmd
  real, intent(in) :: plev_pa(kmd)
  character(len=*), intent(in) :: vname
  real(kind=4), intent(inout) :: field(is:ie,js:je,kmd)
  logical, intent(in) :: valid(is:ie,js:je,kmd)
  integer :: i,j,k,m,kup,kdn
  real :: w,lp

  do j=js,je
    do i=is,ie
      if (.not.any(valid(i,j,1:kmd))) then
        write(*,'(a,a,a,3i8)') 'GRIB2 ',trim(vname),' has no valid vertical data at PE,i,j=', &
                               mpp_pe(),i,j
        call mpp_error(FATAL,'GRIB2 bitmap leaves an entirely missing target column')
      endif
      do k=1,kmd
        if (.not.valid(i,j,k)) then
          kup=0; kdn=0
          do m=k-1,1,-1
            if (valid(i,j,m)) then
              kup=m
              exit
            endif
          enddo
          do m=k+1,kmd
            if (valid(i,j,m)) then
              kdn=m
              exit
            endif
          enddo
          if (kup>0 .and. kdn>0) then
            lp=log(plev_pa(k))
            w=(lp-log(plev_pa(kup)))/(log(plev_pa(kdn))-log(plev_pa(kup)))
            field(i,j,k)=(1.0-w)*field(i,j,kup)+w*field(i,j,kdn)
          else if (kup>0) then
            ! Typical below-ground bitmap case: use the nearest valid level above.
            field(i,j,k)=field(i,j,kup)
          else
            field(i,j,k)=field(i,j,kdn)
          endif
        endif
      enddo
    enddo
  enddo
end subroutine fill_grib2_vertical_missing

subroutine grib2_read_field_open(lugb,category,param,plev_pa,field,valid_mask,lon_deg,lat_deg,return_grid,new_file,required)
  use grib_mod, only: gribfield
  use mpp_mod,  only: mpp_error, FATAL
  implicit none
  integer, intent(in) :: lugb,category,param,plev_pa
  logical, intent(in) :: return_grid,new_file,required
  real, allocatable, intent(inout) :: field(:,:),valid_mask(:,:),lon_deg(:),lat_deg(:)

  type(gribfield) :: gfld
  integer :: jdisc,jpdtn,jgdtn,jskp,kfound,iret
  integer :: jids(200),jpdt(200),jgdt(200)
  integer :: nx_raw,ny,nx_out,scan,i,j,idx,iw,nfirst
  integer :: sf,sv
  real :: unit_deg,lon1,lat1,lon2,lat2,di,dj,si,sj,span,pfound,tol_pa,tol_lon
  real, allocatable :: raw(:,:),work(:,:),work_valid(:,:),tmp(:,:),tmp_valid(:,:),lon_work(:)
  logical :: ineg,jpos,jconsecutive,alternate,level_found

  if (new_file) grib2_time_match_set=.false.

  jdisc = 0
  jpdtn = 0
  jgdtn = -1
  jids  = -9999
  jpdt  = -9999
  jgdt  = -9999
  if (grib2_time_match_set) then
    jids(5:11)=grib2_time_jids(:)
    jpdt(8)=grib2_time_unit
    jpdt(9)=grib2_forecast_time
  endif
  jpdt(1)  = category
  jpdt(2)  = param
  jpdt(10) = 100

  ! Fast/common case: isobaric level encoded directly in Pa with scale factor 0.
  jpdt(11) = 0
  jpdt(12) = plev_pa
  jskp = 0
  call getgb2(lugb,0,jskp,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt, &
              .true.,kfound,gfld,iret)

  ! Generic fallback: accept any fixed-surface scale factor, then compare the
  ! decoded physical pressure and continue searching until the requested level.
  if (iret /= 0) then
    if (iret /= 99) then
      call mpp_error(FATAL,'grib2_read_field_open: error reading GRIB2/index')
    endif
    jpdt(11)=-9999
    jpdt(12)=-9999
    jskp=0
    level_found=.false.
    do
      call getgb2(lugb,0,jskp,jdisc,jids,jpdtn,jpdt,jgdtn,jgdt, &
                  .true.,kfound,gfld,iret)
      if (iret /= 0) exit
      if (gfld%ipdtlen >= 12) then
        sf=gfld%ipdtmpl(11)
        sv=gfld%ipdtmpl(12)
        if (sf > -100 .and. sf < 100) then
          pfound=real(sv)*10.0**(-real(sf))
          tol_pa=max(0.5,1.e-6*real(plev_pa))
          if (abs(pfound-real(plev_pa)) <= tol_pa) then
            level_found=.true.
            exit
          endif
        endif
      endif
      jskp=kfound
      call gf_free(gfld)
    enddo
    if (.not.level_found) then
      if (.not.required .and. iret == 99) then
        if (allocated(field)) field=0.0
        if (allocated(valid_mask)) valid_mask=0.0
        return
      endif
      write(*,'(a,4i10)') 'GRIB2 field not found: cat,param,type,p(Pa)=', &
                           category,param,100,plev_pa
      call mpp_error(FATAL,'grib2_read_field_open: getgb2 failed')
    endif
  endif

  if (.not.grib2_time_match_set) then
    if (gfld%idsectlen >= 11 .and. gfld%ipdtlen >= 9) then
      grib2_time_jids(:)=gfld%idsect(5:11)
      grib2_time_unit=gfld%ipdtmpl(8)
      grib2_forecast_time=gfld%ipdtmpl(9)
      grib2_time_match_set=.true.
    endif
  endif

  if (gfld%igdtnum /= 0) then
    call gf_free(gfld)
    call mpp_error(FATAL,'GRIB2 reader currently requires regular lat-lon grid template 3.0')
  endif
  if (gfld%igdtlen < 19) then
    call gf_free(gfld)
    call mpp_error(FATAL,'GRIB2 grid template is too short')
  endif

  if (return_grid) then
    grib2_grid_signature(1:19)=gfld%igdtmpl(1:19)
    grib2_grid_signature_set=.true.
  else
    if (.not.grib2_grid_signature_set) then
      call gf_free(gfld)
      call mpp_error(FATAL,'GRIB2 source-grid signature was not initialized')
    endif
    if (any(gfld%igdtmpl(8:19) /= grib2_grid_signature(8:19))) then
      call gf_free(gfld)
      call mpp_error(FATAL,'GRIB2 grid definition changed between variables/levels')
    endif
  endif

  nx_raw = gfld%igdtmpl(8)
  ny = gfld%igdtmpl(9)
  if (gfld%ngrdpts /= nx_raw*ny) then
    call gf_free(gfld)
    call mpp_error(FATAL,'GRIB2 grid is not a regular rectangular grid')
  endif

  scan = gfld%igdtmpl(19)
  ineg         = btest(scan,7)
  jpos         = btest(scan,6)
  jconsecutive = btest(scan,5)
  alternate    = btest(scan,4)
  if (jconsecutive .or. alternate) then
    call gf_free(gfld)
    call mpp_error(FATAL,'GRIB2 reader supports only i-consecutive, non-alternating scan order')
  endif

  allocate(raw(nx_raw,ny),work(nx_raw,ny),work_valid(nx_raw,ny))
  raw = reshape(gfld%fld(1:nx_raw*ny),(/nx_raw,ny/))
  work_valid=1.0
  if (gfld%ibmap /= 255 .and. associated(gfld%bmap)) then
    do j=1,ny
      do i=1,nx_raw
        idx=(j-1)*nx_raw+i
        if (.not.gfld%bmap(idx)) work_valid(i,j)=0.0
      enddo
    enddo
  endif

  ! Canonicalize the field and bitmap to increasing i and increasing j.
  if (ineg) then
    do j=1,ny
      do i=1,nx_raw
        work(i,j)=raw(nx_raw+1-i,j)
      enddo
    enddo
    allocate(tmp_valid(nx_raw,ny))
    tmp_valid=work_valid(nx_raw:1:-1,:)
    work_valid=tmp_valid
    deallocate(tmp_valid)
  else
    work=raw
  endif
  if (.not.jpos) then
    allocate(tmp(nx_raw,ny),tmp_valid(nx_raw,ny))
    do j=1,ny
      tmp(:,j)=work(:,ny+1-j)
      tmp_valid(:,j)=work_valid(:,ny+1-j)
    enddo
    work=tmp
    work_valid=tmp_valid
    deallocate(tmp,tmp_valid)
  endif

  if (return_grid) then
    grib2_raw_im=nx_raw
    grib2_raw_jm=ny

    ! For GDT 3.0, basic angle 0 means the standard 1.e-6 degree unit;
    ! subdivisions is normally encoded missing and must not be required to 0.
    if (gfld%igdtmpl(10) == 0) then
      unit_deg=1.0e-6
    else if (gfld%igdtmpl(11) > 0 .and. &
             gfld%igdtmpl(11) < huge(gfld%igdtmpl(11))) then
      unit_deg=real(gfld%igdtmpl(10))/real(gfld%igdtmpl(11))
    else
      call gf_free(gfld)
      call mpp_error(FATAL,'GRIB2 reader found invalid lat-lon angular scaling')
    endif

    lon1=real(gfld%igdtmpl(13))*unit_deg
    lat1=real(gfld%igdtmpl(12))*unit_deg
    lon2=real(gfld%igdtmpl(16))*unit_deg
    lat2=real(gfld%igdtmpl(15))*unit_deg
    si=merge(-1.0,1.0,ineg)
    sj=merge( 1.0,-1.0,jpos)

    if (gfld%igdtmpl(17) > 0 .and. &
        gfld%igdtmpl(17) < huge(gfld%igdtmpl(17))) then
      di=abs(real(gfld%igdtmpl(17))*unit_deg)
    else
      if (nx_raw <= 1) call mpp_error(FATAL,'GRIB2 cannot derive Di with Ni <= 1')
      if (ineg) then
        span=modulo(lon1-lon2,360.0)
      else
        span=modulo(lon2-lon1,360.0)
      endif
      if (span < 1.e-8) span=360.0
      di=span/real(nx_raw-1)
    endif
    if (gfld%igdtmpl(18) > 0 .and. &
        gfld%igdtmpl(18) < huge(gfld%igdtmpl(18))) then
      dj=abs(real(gfld%igdtmpl(18))*unit_deg)
    else
      if (ny <= 1) call mpp_error(FATAL,'GRIB2 cannot derive Dj with Nj <= 1')
      dj=abs(lat2-lat1)/real(ny-1)
    endif
    if (di <= 0.0 .or. dj <= 0.0) then
      call gf_free(gfld)
      call mpp_error(FATAL,'GRIB2 reader found invalid/zero Di or Dj')
    endif

    if (allocated(lat_deg)) deallocate(lat_deg)
    allocate(lon_work(nx_raw),lat_deg(ny))
    do i=1,nx_raw
      lon_work(i)=modulo(lon1+si*real(i-1)*di,360.0)
    enddo
    if (ineg) lon_work=lon_work(nx_raw:1:-1)
    do j=1,ny
      lat_deg(j)=lat1+sj*real(j-1)*dj
    enddo
    if (.not.jpos) lat_deg=lat_deg(ny:1:-1)

    ! A global grid may contain both 0 and 360 degrees.  After modulo these
    ! are duplicate endpoints; drop the last canonical column before rotation.
    tol_lon=max(1.e-6,1.e-3*di)
    grib2_drop_duplicate_lon=.false.
    if (nx_raw>1) then
      if (abs(di*real(nx_raw-1)-360.0) <= max(tol_lon,1.e-4) .and. &
          abs(modulo(lon_work(nx_raw)-lon_work(1)+180.0,360.0)-180.0) <= tol_lon) then
        grib2_drop_duplicate_lon=.true.
      endif
    endif
    nx_out=nx_raw
    if (grib2_drop_duplicate_lon) nx_out=nx_raw-1

    ! Rotate a wrapped global longitude array once, and save the same transform
    ! for every subsequent variable and pressure level.
    grib2_lon_start=1
    do i=2,nx_out
      if (lon_work(i) <= lon_work(i-1)) then
        grib2_lon_start=i
        exit
      endif
    enddo

    if (allocated(field)) deallocate(field)
    if (allocated(valid_mask)) deallocate(valid_mask)
    if (allocated(lon_deg)) deallocate(lon_deg)
    allocate(field(nx_out,ny),valid_mask(nx_out,ny),lon_deg(nx_out))

    if (grib2_lon_start > 1) then
      iw=grib2_lon_start
      nfirst=nx_out-iw+1
      field(1:nfirst,:)=work(iw:nx_out,:)
      field(nfirst+1:nx_out,:)=work(1:iw-1,:)
      valid_mask(1:nfirst,:)=work_valid(iw:nx_out,:)
      valid_mask(nfirst+1:nx_out,:)=work_valid(1:iw-1,:)
      lon_deg(1:nfirst)=lon_work(iw:nx_out)
      lon_deg(nfirst+1:nx_out)=lon_work(1:iw-1)
    else
      field=work(1:nx_out,:)
      valid_mask=work_valid(1:nx_out,:)
      lon_deg=lon_work(1:nx_out)
    endif

    do i=2,nx_out
      if (lon_deg(i) <= lon_deg(i-1)) then
        call gf_free(gfld)
        call mpp_error(FATAL,'GRIB2 longitude array is not strictly increasing')
      endif
    enddo
    do j=2,ny
      if (lat_deg(j) <= lat_deg(j-1)) then
        call gf_free(gfld)
        call mpp_error(FATAL,'GRIB2 latitude array is not strictly increasing')
      endif
    enddo
    deallocate(lon_work)
  else
    if (nx_raw /= grib2_raw_im .or. ny /= grib2_raw_jm) then
      call gf_free(gfld)
      call mpp_error(FATAL,'GRIB2 source grid changed between files/fields')
    endif
    nx_out=nx_raw
    if (grib2_drop_duplicate_lon) nx_out=nx_raw-1
    if (nx_out /= im .or. ny /= jm) then
      call gf_free(gfld)
      call mpp_error(FATAL,'GRIB2 canonical source grid changed between files/fields')
    endif
    if (.not.allocated(field)) allocate(field(nx_out,ny))
    if (.not.allocated(valid_mask)) allocate(valid_mask(nx_out,ny))
    if (size(field,1)/=nx_out .or. size(field,2)/=ny) then
      deallocate(field); allocate(field(nx_out,ny))
    endif
    if (size(valid_mask,1)/=nx_out .or. size(valid_mask,2)/=ny) then
      deallocate(valid_mask); allocate(valid_mask(nx_out,ny))
    endif
    if (grib2_lon_start > 1) then
      iw=grib2_lon_start
      nfirst=nx_out-iw+1
      field(1:nfirst,:)=work(iw:nx_out,:)
      field(nfirst+1:nx_out,:)=work(1:iw-1,:)
      valid_mask(1:nfirst,:)=work_valid(iw:nx_out,:)
      valid_mask(nfirst+1:nx_out,:)=work_valid(1:iw-1,:)
    else
      field=work(1:nx_out,:)
      valid_mask=work_valid(1:nx_out,:)
    endif
  endif

  deallocate(raw,work,work_valid)
  call gf_free(gfld)
end subroutine grib2_read_field_open

!>@brief Interpolate pressure-level GRIB2 T and q to the actual FV3 full-level pressure.
!! The GRIB2 values are point values on fixed isobaric surfaces, so no
!! artificial source layers are constructed.  The target FV3 interface
!! pressures are built from the actual model delp, as in remap_tq/remap_uv,
!! and the full-level pressure uses the FV3 logarithmic-mean definition.
!! The GRIB2 profile is then linearly interpolated in log(p). Pressure
!! surfaces below the actual FV3 surface are excluded, and extrapolation is
!! bounded to one source-level interval at the top and usable bottom of a column.
 subroutine remap_tq_grib2(npz, ak, bk, ps, delp, t, q, &
                           kmd, plev_pa, ta, qa, zvir, ptop)
  integer, intent(in):: npz, kmd
  real,    intent(in):: zvir, ptop
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(kmd):: plev_pa
  real,    intent(inout), dimension(is:ie,js:je):: ps
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4), intent(in), dimension(is:ie,js:je,kmd):: ta, qa
  real(KIND=4), intent(out), dimension(is:ie,js:je,npz):: t
! on  input "ta" is virtual temperature
! on output "t" is virtual temperature
  real(KIND=4), intent(out), dimension(is:ie,js:je,npz):: q
! local:
  real, dimension(kmd):: lplev
  real, dimension(is:ie,npz+1):: pe1
  real, dimension(is:ie,npz):: pfull
  real :: lpfull, w
  integer i,j,k,m,kbot_src
  real, parameter :: pe_top_eps = 1.e-9

  if (kmd < 2) then
     call mpp_error(FATAL,'remap_tq_grib2: at least two pressure levels are required')
  endif
  do k=1,kmd
     if (plev_pa(k) <= 0.0) then
        call mpp_error(FATAL,'remap_tq_grib2: source pressure levels must be positive')
     endif
     lplev(k) = log(plev_pa(k))
  enddo
  do k=2,kmd
     if (plev_pa(k) <= plev_pa(k-1)) then
        call mpp_error(FATAL,'remap_tq_grib2: source pressure levels must increase downward')
     endif
  enddo

  if ( .not. nudge_q .and. .not. climate_nudging ) q = 0.0

  do 5000 j=js,je

!------
! Model
!------
! Build the actual FV3 interface pressure from the current delp, exactly as
! in the existing remap_tq/remap_uv target-grid construction.
     do i=is,ie
        pe1(i,1) = max(pe_top_eps,ak(1))
     enddo
     do k=1,npz
        do i=is,ie
           pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
           if (pe1(i,k+1) <= pe1(i,k)) then
              call mpp_error(FATAL,'remap_tq_grib2: non-positive model layer pressure thickness')
           endif
! FV3 full-level pressure: logarithmic mean of the two interface pressures.
           pfull(i,k) = (pe1(i,k+1)-pe1(i,k)) / &
                        log(pe1(i,k+1)/pe1(i,k))
        enddo
     enddo
     do i=is,ie
        ps(i,j) = pe1(i,npz+1)
     enddo

!------
! map q
!------
   if ( nudge_q .or. climate_nudging ) then
      do i=is,ie
         ! Do not use pressure surfaces that are below the actual FV3 surface.
         ! This cannot diagnose source-model underground values when the GRIB2
         ! file has no source surface pressure, but it prevents target-terrain
         ! interpolation from using levels that are underground in FV3.
         kbot_src=1
         do m=2,kmd
            if (plev_pa(m) <= pe1(i,npz+1)) kbot_src=m
         enddo
         kbot_src=max(2,kbot_src)
         do k=1,npz
            lpfull = log(pfull(i,k))
            if (pfull(i,k) <= plev_pa(1)) then
               w = max(-1.0,(lpfull-lplev(1))/(lplev(2)-lplev(1)))
               q(i,j,k) = (1.0-w)*qa(i,j,1) + w*qa(i,j,2)
            else if (pfull(i,k) >= plev_pa(kbot_src)) then
               w = min(2.0,(lpfull-lplev(kbot_src-1))/ &
                           (lplev(kbot_src)-lplev(kbot_src-1)))
               q(i,j,k) = (1.0-w)*qa(i,j,kbot_src-1) + w*qa(i,j,kbot_src)
            else
               do m=1,kbot_src-1
                  if (pfull(i,k) >= plev_pa(m) .and. &
                      pfull(i,k) <= plev_pa(m+1)) then
                     w = (lpfull-lplev(m))/(lplev(m+1)-lplev(m))
                     q(i,j,k) = (1.0-w)*qa(i,j,m) + w*qa(i,j,m+1)
                     exit
                  endif
               enddo
            endif
            q(i,j,k)=max(0.0,q(i,j,k))
         enddo
      enddo
   endif

!------
! map t
!------
   do i=is,ie
      kbot_src=1
      do m=2,kmd
         if (plev_pa(m) <= pe1(i,npz+1)) kbot_src=m
      enddo
      kbot_src=max(2,kbot_src)
      do k=1,npz
         lpfull = log(pfull(i,k))
         if (pfull(i,k) <= plev_pa(1)) then
            w = max(-1.0,(lpfull-lplev(1))/(lplev(2)-lplev(1)))
            t(i,j,k) = (1.0-w)*ta(i,j,1) + w*ta(i,j,2)
         else if (pfull(i,k) >= plev_pa(kbot_src)) then
            w = min(2.0,(lpfull-lplev(kbot_src-1))/ &
                        (lplev(kbot_src)-lplev(kbot_src-1)))
            t(i,j,k) = (1.0-w)*ta(i,j,kbot_src-1) + w*ta(i,j,kbot_src)
         else
            do m=1,kbot_src-1
               if (pfull(i,k) >= plev_pa(m) .and. &
                   pfull(i,k) <= plev_pa(m+1)) then
                  w = (lpfull-lplev(m))/(lplev(m+1)-lplev(m))
                  t(i,j,k) = (1.0-w)*ta(i,j,m) + w*ta(i,j,m+1)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

5000 continue

 end subroutine remap_tq_grib2

!>@brief Interpolate pressure-level GRIB2 u and v to the actual FV3 full-level pressure.
!! Horizontal interpolation has already placed each isobaric GRIB2 field on
!! the local FV3 A-grid.  This routine builds the actual FV3 pressure column
!! from delp, evaluates the FV3 logarithmic-mean full-level pressure, excludes
!! pressure surfaces below the FV3 surface, and uses bounded interpolation/
!! extrapolation of the pressure-level wind profile in log(p).
 subroutine remap_uv_grib2(npz, ak, bk, ps, delp, u, v, &
                           kmd, plev_pa, u0, v0, ptop)
  integer, intent(in):: npz, kmd
  real,    intent(in):: ptop
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(kmd):: plev_pa
  real,    intent(inout):: ps(is:ie,js:je)
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4), intent(in), dimension(is:ie,js:je,kmd):: u0, v0
  real(KIND=4), intent(inout), dimension(is:ie,js:je,npz):: u, v
! local:
  real, dimension(kmd):: lplev
  real, dimension(is:ie,npz+1):: pe1
  real, dimension(is:ie,npz):: pfull
  real :: lpfull, w
  integer i,j,k,m,kbot_src
  real, parameter :: pe_top_eps = 1.e-9

  if (kmd < 2) then
     call mpp_error(FATAL,'remap_uv_grib2: at least two pressure levels are required')
  endif
  do k=1,kmd
     if (plev_pa(k) <= 0.0) then
        call mpp_error(FATAL,'remap_uv_grib2: source pressure levels must be positive')
     endif
     lplev(k) = log(plev_pa(k))
  enddo
  do k=2,kmd
     if (plev_pa(k) <= plev_pa(k-1)) then
        call mpp_error(FATAL,'remap_uv_grib2: source pressure levels must increase downward')
     endif
  enddo

  do 5000 j=js,je
!------
! Model
!------
! Use the same actual target interface-pressure construction as remap_uv.
     do i=is,ie
        pe1(i,1) = max(pe_top_eps,ak(1))
     enddo
     do k=1,npz
        do i=is,ie
           pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
           if (pe1(i,k+1) <= pe1(i,k)) then
              call mpp_error(FATAL,'remap_uv_grib2: non-positive model layer pressure thickness')
           endif
! FV3 full-level pressure: logarithmic mean of the two interface pressures.
           pfull(i,k) = (pe1(i,k+1)-pe1(i,k)) / &
                        log(pe1(i,k+1)/pe1(i,k))
        enddo
     enddo
     do i=is,ie
        ps(i,j) = pe1(i,npz+1)
     enddo
!------
! map u
!------
     do i=is,ie
        kbot_src=1
        do m=2,kmd
           if (plev_pa(m) <= pe1(i,npz+1)) kbot_src=m
        enddo
        kbot_src=max(2,kbot_src)
        do k=1,npz
           lpfull = log(pfull(i,k))
           if (pfull(i,k) <= plev_pa(1)) then
              w = max(-1.0,(lpfull-lplev(1))/(lplev(2)-lplev(1)))
              u(i,j,k) = (1.0-w)*u0(i,j,1) + w*u0(i,j,2)
           else if (pfull(i,k) >= plev_pa(kbot_src)) then
              w = min(2.0,(lpfull-lplev(kbot_src-1))/ &
                          (lplev(kbot_src)-lplev(kbot_src-1)))
              u(i,j,k) = (1.0-w)*u0(i,j,kbot_src-1) + w*u0(i,j,kbot_src)
           else
              do m=1,kbot_src-1
                 if (pfull(i,k) >= plev_pa(m) .and. &
                     pfull(i,k) <= plev_pa(m+1)) then
                    w = (lpfull-lplev(m))/(lplev(m+1)-lplev(m))
                    u(i,j,k) = (1.0-w)*u0(i,j,m) + w*u0(i,j,m+1)
                    exit
                 endif
              enddo
           endif
        enddo
     enddo
!------
! map v
!------
     do i=is,ie
        kbot_src=1
        do m=2,kmd
           if (plev_pa(m) <= pe1(i,npz+1)) kbot_src=m
        enddo
        kbot_src=max(2,kbot_src)
        do k=1,npz
           lpfull = log(pfull(i,k))
           if (pfull(i,k) <= plev_pa(1)) then
              w = max(-1.0,(lpfull-lplev(1))/(lplev(2)-lplev(1)))
              v(i,j,k) = (1.0-w)*v0(i,j,1) + w*v0(i,j,2)
           else if (pfull(i,k) >= plev_pa(kbot_src)) then
              w = min(2.0,(lpfull-lplev(kbot_src-1))/ &
                          (lplev(kbot_src)-lplev(kbot_src-1)))
              v(i,j,k) = (1.0-w)*v0(i,j,kbot_src-1) + w*v0(i,j,kbot_src)
           else
              do m=1,kbot_src-1
                 if (pfull(i,k) >= plev_pa(m) .and. &
                     pfull(i,k) <= plev_pa(m+1)) then
                    w = (lpfull-lplev(m))/(lplev(m+1)-lplev(m))
                    v(i,j,k) = (1.0-w)*v0(i,j,m) + w*v0(i,j,m+1)
                    exit
                 endif
              enddo
           endif
        enddo
     enddo
5000 continue

 end subroutine remap_uv_grib2

 subroutine remap_coef(agrid)

  real, intent(IN), dimension(isd:ied,jsd:jed,2) :: agrid

! local:
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1
  integer i,j, i1, i2, jc, i0, j0

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to cubed sphere cell center
  do 5000 j=js,je

     do i=is,ie

       if ( agrid(i,j,1)>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)-lon(im)) * rdlon(im)
       elseif ( agrid(i,j,1)<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (agrid(i,j,1)+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=1,im-1
            if ( agrid(i,j,1)>=lon(i0) .and. agrid(i,j,1)<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               a1 = (agrid(i,j,1)-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue

       if ( agrid(i,j,2)<lat(1) ) then
            jc = 1
            b1 = 0.
       elseif ( agrid(i,j,2)>lat(jm) ) then
            jc = jm-1
            b1 = 1.
       else
          do j0=1,jm-1
          if ( agrid(i,j,2)>=lat(j0) .and. agrid(i,j,2)<=lat(j0+1) ) then
               jc = j0
               b1 = (agrid(i,j,2)-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
       endif
222    continue

       if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
            write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
       endif

       s2c(i,j,1) = (1.-a1) * (1.-b1)
       s2c(i,j,2) =     a1  * (1.-b1)
       s2c(i,j,3) =     a1  *     b1
       s2c(i,j,4) = (1.-a1) *     b1
       id1(i,j) = i1
       id2(i,j) = i2
       jdc(i,j) = jc
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine remap_coef

#ifndef DYCORE_SOLO
 subroutine ncep2fms( sst )
  real(kind=4), intent(in):: sst(im,jm)
! local:
  real :: rdlon(im)
  real :: rdlat(jm)
  real:: a1, b1
  real:: delx, dely
  real:: xc, yc    ! "data" location
  real:: c1, c2, c3, c4
  integer i,j, i1, i2, jc, i0, j0, it, jt

  do i=1,im-1
     rdlon(i) = 1. / (lon(i+1) - lon(i))
  enddo
     rdlon(im) = 1. / (lon(1) + 2.*pi - lon(im))

  do j=1,jm-1
     rdlat(j) = 1. / (lat(j+1) - lat(j))
  enddo

! * Interpolate to "FMS" 1x1 SST data grid
! lon: 0.5, 1.5, ..., 359.5
! lat: -89.5, -88.5, ... , 88.5, 89.5

  delx = 360./real(i_sst)
  dely = 180./real(j_sst)

  jt = 1
  do 5000 j=1,j_sst

     yc = (-90. + dely * (0.5+real(j-1)))  * deg2rad
     if ( yc<lat(1) ) then
            jc = 1
            b1 = 0.
     elseif ( yc>lat(jm) ) then
            jc = jm-1
            b1 = 1.
     else
          do j0=jt,jm-1
          if ( yc>=lat(j0) .and. yc<=lat(j0+1) ) then
               jc = j0
               jt = j0
               b1 = (yc-lat(jc)) * rdlat(jc)
               go to 222
          endif
          enddo
     endif
222  continue
     it = 1

     do i=1,i_sst
        xc = delx * (0.5+real(i-1)) * deg2rad
       if ( xc>lon(im) ) then
            i1 = im;     i2 = 1
            a1 = (xc-lon(im)) * rdlon(im)
       elseif ( xc<lon(1) ) then
            i1 = im;     i2 = 1
            a1 = (xc+2.*pi-lon(im)) * rdlon(im)
       else
            do i0=it,im-1
            if ( xc>=lon(i0) .and. xc<=lon(i0+1) ) then
               i1 = i0;  i2 = i0+1
               it = i0
               a1 = (xc-lon(i1)) * rdlon(i0)
               go to 111
            endif
            enddo
       endif
111    continue

!      if ( a1<0.0 .or. a1>1.0 .or.  b1<0.0 .or. b1>1.0 ) then
!           write(*,*) 'gid=', mpp_pe(), i,j,a1, b1
!      endif
       c1 = (1.-a1) * (1.-b1)
       c2 =     a1  * (1.-b1)
       c3 =     a1  *     b1
       c4 = (1.-a1) *     b1
! Interpolated surface pressure
       sst_ncep(i,j) = c1*sst(i1,jc  ) + c2*sst(i2,jc  ) +    &
                       c3*sst(i2,jc+1) + c4*sst(i1,jc+1)
     enddo   !i-loop
5000 continue   ! j-loop

 end subroutine ncep2fms
#endif

 subroutine get_int_hght(npz, ak, bk, ps, delp, ps0, tv)
  integer, intent(in):: npz
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in), dimension(is:ie,js:je):: ps, ps0
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4),  intent(in), dimension(is:ie,js:je,km):: tv    ! virtual temperature
! local:
  real, dimension(is:ie,km+1):: pn0
  integer i,j,k

  ! NOTE: this routine is only called when nudge_ps=.true.  For nudge_option=1
  ! that path also needs gz3 allocated (done conditionally in fv_nwp_nudge_init)
  ! and ak0/bk0 to reflect the analysis vertical grid.  Since the regional
  ! analysis stores delp directly, ak0 is set as a sentinel (see init); if
  ! nudge_ps is ever enabled for nudge_option=1, ak0/bk0 must be loaded from
  ! the analysis file's hybrid coefficients before this routine will be correct.
!$OMP parallel do default(none) shared(is,ie,js,je,km,ak0,bk0,ps0,gz3,gz0,tv) &
!$OMP                          private(pn0)
  do 5000 j=js,je

     do k=1,km+1
        do i=is,ie
           pn0(i,k) = log( ak0(k) + bk0(k)*ps0(i,j) )
        enddo
     enddo
     do i=is,ie
        gz3(i,j,km+1) = gz0(i,j)   ! Data Surface geopotential
     enddo

     do k=km,1,-1
        do i=is,ie
           gz3(i,j,k) = gz3(i,j,k+1) + rdgas*tv(i,j,k)*(pn0(i,k+1)-pn0(i,k))
        enddo
     enddo

5000 continue

 end subroutine get_int_hght

 subroutine remap_tq( npz, ak,  bk,  ps, delp,  t,  q,  &
                      kmd, ps0, ta, qa, zvir, ptop)
  integer, intent(in):: npz, kmd
  real,    intent(in):: zvir, ptop
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(in),    dimension(is:ie,js:je):: ps0
  real,    intent(inout), dimension(is:ie,js:je):: ps
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4),    intent(in), dimension(is:ie,js:je,kmd):: ta, qa
  real(KIND=4),    intent(out), dimension(is:ie,js:je,npz):: t
! on  input "ta" is virtual temperature
! on output "t" is virtual temperature
  real(KIND=4),    intent(out), dimension(is:ie,js:je,npz):: q
! local:
  real, dimension(is:ie,kmd):: tp, qp
  real, dimension(is:ie,kmd+1):: pe0, pn0
  real, dimension(is:ie,npz):: qn1
  real, dimension(is:ie,npz+1):: pe1, pn1
  integer i,j,k

  do 5000 j=js,je

     do k=1,kmd+1
        do i=is,ie
           pe0(i,k) = ak0(k) + bk0(k)*ps0(i,j)
           pn0(i,k) = log(pe0(i,k))
       enddo
     enddo
!------
! Model
!------
     do i=is,ie
        pe1(i,1) = ak(1)
     enddo
     do k=1,npz
       do i=is,ie
          pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
       enddo
     enddo
     do i=is,ie
        ps(i,j) = pe1(i,npz+1)
     enddo
     do k=1,npz+1
        do i=is,ie
           pn1(i,k) = log(pe1(i,k))
        enddo
     enddo

   if ( nudge_q ) then
        do k=1,kmd
           do i=is,ie
              qp(i,k) = qa(i,j,k)
           enddo
        enddo
        call mappm(kmd, pe0, qp, npz, pe1, qn1, is,ie, 0, kord_data, ptop)
        do k=1,npz
           do i=is,ie
              q(i,j,k) = qn1(i,k)
           enddo
        enddo
   endif

   do k=1,kmd
      do i=is,ie
         tp(i,k) = ta(i,j,k)
      enddo
   enddo
   call mappm(kmd, pn0, tp, npz, pn1, qn1, is,ie, 1, kord_data, ptop)

   do k=1,npz
      do i=is,ie
         t(i,j,k) = qn1(i,k)
      enddo
   enddo

5000 continue

 end subroutine remap_tq

 subroutine remap_uv(npz, ak, bk, ps, delp, u, v, kmd, ps0, u0, v0, ptop)
  integer, intent(in):: npz
  real,    intent(IN):: ptop
  real,    intent(in):: ak(npz+1), bk(npz+1)
  real,    intent(inout):: ps(is:ie,js:je)
  real, intent(in), dimension(isd:ied,jsd:jed,npz):: delp
  real(KIND=4),    intent(inout), dimension(is:ie,js:je,npz):: u, v
!
  integer, intent(in):: kmd
  real,    intent(in):: ps0(is:ie,js:je)
  real(KIND=4),    intent(in), dimension(is:ie,js:je,kmd):: u0, v0
!
! local:
  real, dimension(is:ie,kmd+1):: pe0
  real, dimension(is:ie,npz+1):: pe1
  real, dimension(is:ie,kmd):: qt
  real, dimension(is:ie,npz):: qn1
  integer i,j,k

  do 5000 j=js,je
!------
! Data
!------
     do k=1,kmd+1
       do i=is,ie
          pe0(i,k) = ak0(k) + bk0(k)*ps0(i,j)
       enddo
     enddo
!------
! Model
!------
     do i=is,ie
        pe1(i,1) = ak(1)
     enddo
     do k=1,npz
       do i=is,ie
          pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
       enddo
     enddo
     do i=is,ie
        ps(i,j) = pe1(i,npz+1)
     enddo
!------
! map u
!------
      do k=1,kmd
         do i=is,ie
            qt(i,k) = u0(i,j,k)
         enddo
      enddo
      call mappm(kmd, pe0, qt, npz, pe1, qn1, is,ie, -1, kord_data, ptop)
      do k=1,npz
         do i=is,ie
            u(i,j,k) = qn1(i,k)
         enddo
      enddo
!------
! map v
!------
      do k=1,kmd
         do i=is,ie
            qt(i,k) = v0(i,j,k)
         enddo
      enddo
      call mappm(kmd, pe0, qt, npz, pe1, qn1, is,ie, -1, kord_data, ptop)
      do k=1,npz
         do i=is,ie
            v(i,j,k) = qn1(i,k)
         enddo
      enddo
5000 continue

 end subroutine remap_uv

!>@brief Regional remap of T and q using dp0 (analysis delp) for source grid.
!! Identical to remap_tq except pe0 is built by cumulative sum of dp0
!! rather than ak0+bk0*ps0, avoiding the ak0 size/content mismatch that
!! occurs when kmd (analysis levels) differs from npz (model levels).
 subroutine remap_tq_regional(npz, ak, bk, ps, delp, t, q, &
                               kmd, dp0, ta, qa, zvir, ptop)
  integer,      intent(in)    :: npz, kmd
  real,         intent(in)    :: zvir, ptop
  real,         intent(in)    :: ak(npz+1), bk(npz+1)
  real,         intent(inout), dimension(is:ie,js:je)         :: ps
  real,         intent(in),    dimension(isd:ied,jsd:jed,npz) :: delp
  real(KIND=4), intent(in),    dimension(is:ie,js:je,kmd)     :: dp0, ta, qa
  real(KIND=4), intent(out),   dimension(is:ie,js:je,npz)     :: t, q
! local:
  real, dimension(is:ie,kmd)   :: tp, qp
  real, dimension(is:ie,kmd+1) :: pe0, pn0
  real, dimension(is:ie,npz)   :: qn1
  real, dimension(is:ie,npz+1) :: pe1, pn1
  integer :: i, j, k
  real, parameter :: pe_top_eps = 1.e-9

  if ( .not. nudge_q .and. .not. climate_nudging ) q = 0.0

  do 5000 j = js, je
!------
! Analysis (source) grid: pe0(i,1)=1e-9 Pa approximating the true vacuum
! top (ak0(1)=0). This is negligibly small vs the model top ak(1)~200 Pa
! and avoids log(0) for pn0 without any special clamping.
! Cumulative sum gives pe0(i,km+1)=ps0 (correct analysis surface pressure).
!------
     do i = is, ie
        pe0(i,1) = pe_top_eps
        pn0(i,1) = log(pe0(i,1))
     enddo
     do k = 1, kmd
        do i = is, ie
           pe0(i,k+1) = pe0(i,k) + real(dp0(i,j,k))
           pn0(i,k+1) = log(pe0(i,k+1))
        enddo
     enddo
!------
! Model (target) grid
!------
     do i = is, ie
        pe1(i,1) = max(pe_top_eps, ak(1))
     enddo
     do k = 1, npz
        do i = is, ie
           pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
        enddo
     enddo
     do i = is, ie
        ps(i,j) = pe1(i,npz+1)
     enddo
     do k = 1, npz+1
        do i = is, ie
           pn1(i,k) = log(pe1(i,k))
        enddo
     enddo

     if ( nudge_q .or. climate_nudging ) then
        do k = 1, kmd
           do i = is, ie
              qp(i,k) = qa(i,j,k)
           enddo
        enddo
        call mappm(kmd, pe0, qp, npz, pe1, qn1, is, ie, 0, kord_data, ptop)
        do k = 1, npz
           do i = is, ie
              q(i,j,k) = qn1(i,k)
           enddo
        enddo
     endif

     do k = 1, kmd
        do i = is, ie
           tp(i,k) = ta(i,j,k)
        enddo
     enddo
     call mappm(kmd, pn0, tp, npz, pn1, qn1, is, ie, 1, kord_data, ptop)
     do k = 1, npz
        do i = is, ie
           t(i,j,k) = qn1(i,k)
        enddo
     enddo

5000 continue

 end subroutine remap_tq_regional

!>@brief Regional remap of u and v using dp0 (analysis delp) for source grid.
!! Identical to remap_uv except pe0 is built by cumulative sum of dp0.
 subroutine remap_uv_regional(npz, ak, bk, ps, delp, u, v, &
                               kmd, dp0, u0, v0, ptop)
  integer,      intent(in)    :: npz, kmd
  real,         intent(in)    :: ptop
  real,         intent(in)    :: ak(npz+1), bk(npz+1)
  real,         intent(inout) :: ps(is:ie,js:je)
  real,         intent(in),    dimension(isd:ied,jsd:jed,npz) :: delp
  real(KIND=4), intent(in),    dimension(is:ie,js:je,kmd)     :: dp0, u0, v0
  real(KIND=4), intent(inout), dimension(is:ie,js:je,npz)     :: u, v
! local:
  real, dimension(is:ie,kmd+1) :: pe0
  real, dimension(is:ie,npz+1) :: pe1
  real, dimension(is:ie,kmd)   :: qt
  real, dimension(is:ie,npz)   :: qn1
  integer :: i, j, k
  real, parameter :: pe_top_eps = 1.e-9

  do 5000 j = js, je
!------
! Analysis (source) grid: pe0(i,1)=1e-9 Pa approximating the true vacuum
! top (ak0(1)=0). Consistent with remap_tq_regional.
! Cumulative sum gives pe0(i,km+1)=ps0 (correct analysis surface pressure).
!------
     do i = is, ie
        pe0(i,1) = pe_top_eps
     enddo
     do k = 1, kmd
        do i = is, ie
           pe0(i,k+1) = pe0(i,k) + real(dp0(i,j,k))
        enddo
     enddo
!------
! Model (target) grid
!------
     do i = is, ie
        pe1(i,1) = max(pe_top_eps, ak(1))
     enddo
     do k = 1, npz
        do i = is, ie
           pe1(i,k+1) = pe1(i,k) + delp(i,j,k)
        enddo
     enddo
     do i = is, ie
        ps(i,j) = pe1(i,npz+1)
     enddo
!------
! map u
!------
     do k = 1, kmd
        do i = is, ie
           qt(i,k) = u0(i,j,k)
        enddo
     enddo
     call mappm(kmd, pe0, qt, npz, pe1, qn1, is, ie, -1, kord_data, ptop)
     do k = 1, npz
        do i = is, ie
           u(i,j,k) = qn1(i,k)
        enddo
     enddo
!------
! map v
!------
     do k = 1, kmd
        do i = is, ie
           qt(i,k) = v0(i,j,k)
        enddo
     enddo
     call mappm(kmd, pe0, qt, npz, pe1, qn1, is, ie, -1, kord_data, ptop)
     do k = 1, npz
        do i = is, ie
           v(i,j,k) = qn1(i,k)
        enddo
     enddo

5000 continue

 end subroutine remap_uv_regional

!>@brief The subroutine 'fv_nwp_nudge_end' terminates the nudging module.
 subroutine fv_nwp_nudge_end
    if ( nudge_option.eq.0 ) then
       deallocate ( ps_dat )
       deallocate (  t_dat )
       deallocate (  q_dat )

       if ( nudge_winds ) then
            deallocate ( u_dat )
            deallocate ( v_dat )
       endif

       deallocate ( s2c )
       deallocate ( id1 )
       deallocate ( id2 )
       deallocate ( jdc )

       deallocate ( ak0 )
       deallocate ( bk0 )
       deallocate ( lat )
       deallocate ( lon )

       deallocate ( gz3 )
       deallocate ( gz0 )
    else
       if (allocated(ps_dat)) deallocate(ps_dat)
       if (allocated(t_dat))  deallocate(t_dat)
       if (allocated(q_dat))  deallocate(q_dat)
       if (allocated(u_dat))  deallocate(u_dat)
       if (allocated(v_dat))  deallocate(v_dat)
       if (allocated(dp_dat)) deallocate(dp_dat)

       if (allocated(s2c)) deallocate(s2c)
       if (allocated(id1)) deallocate(id1)
       if (allocated(id2)) deallocate(id2)
       if (allocated(jdc)) deallocate(jdc)
       if (allocated(gz3)) deallocate(gz3)
       if (allocated(gz0)) deallocate(gz0)
       if (allocated(lat)) deallocate(lat)
       if (allocated(lon)) deallocate(lon)
       if (allocated(ak0)) deallocate(ak0)
       if (allocated(bk0)) deallocate(bk0)

       if (allocated(lat_weight_t)) deallocate(lat_weight_t)
       if (allocated(lat_weight_uv)) deallocate(lat_weight_uv)
       lat_taper_initialized=.false.
       if (allocated(dt_stored)) deallocate(dt_stored)
       if (allocated(du_stored)) deallocate(du_stored)
       if (allocated(dv_stored)) deallocate(dv_stored)
       if (allocated(klist)) deallocate(klist)
       skip_counter=0

       if (allocated(alpha_x_glob)) deallocate(alpha_x_glob)
       if (allocated(alpha_y_glob)) deallocate(alpha_y_glob)
       if (allocated(alpha_x_cube)) deallocate(alpha_x_cube)
       if (allocated(alpha_y_cube)) deallocate(alpha_y_cube)
       rf_alpha_stored=.false.
       call rf_clear_cube_paths()

       grib2_raw_im=0
       grib2_raw_jm=0
       grib2_lon_start=1
       grib2_drop_duplicate_lon=.false.
       grib2_grid_signature_set=.false.
       grib2_time_match_set=.false.
       no_obs=.false.
       analysis_time=.false.
       time_nudge=0.0
       module_is_initialized=.false.
    endif
 end subroutine fv_nwp_nudge_end

 subroutine get_tc_mask(time, mask, agrid)
      real :: slp_mask = 100900.    !< crtical SLP to apply mask
! Input
      type(time_type), intent(in):: time
      real, intent(inout):: mask(is:ie,js:je)
      real(kind=R_GRID), intent(IN), dimension(isd:ied,jsd:jed,2) :: agrid
! local
      real(kind=R_GRID):: pos(2)
      real:: slp_o         !< sea-level pressure (Pa)
      real:: w10_o         !< 10-m wind
      real:: r_vor, p_vor
      real:: dist
      integer n, i, j

    do 5000 n=1,nstorms      ! loop through all storms
!----------------------------------------
! Obtain slp observation
!----------------------------------------
      call get_slp_obs(time, nobs_tc(n), x_obs(1,n), y_obs(1,n), wind_obs(1,n),  mslp_obs(1,n), mslp_out(1,n), rad_out(1,n),   &
                       time_tc(1,n), pos(1), pos(2), w10_o, slp_o, r_vor, p_vor)

      if ( slp_o<880.E2 .or. slp_o>min(slp_env,slp_mask) .or. abs(pos(2))*rad2deg>40. ) goto 5000  ! next storm

      if ( r_vor < 30.E3 ) then
           r_vor = r_min + (slp_env-slp_o)/20.E2*r_inc   ! radius of influence
      endif

      do j=js, je
         do i=is, ie
            dist = great_circle_dist(pos, agrid(i,j,1:2), radius)
            if( dist < 6.*r_vor  ) then
                mask(i,j) = mask(i,j) * ( 1. - mask_fac*exp(-(0.5*dist/r_vor)**2)*min(1.,(slp_env-slp_o)/10.E2) )
            endif
         enddo             ! i-loop
      enddo                ! end j-loop

5000 continue

 end subroutine get_tc_mask

!>@brief The subroutine 'breed_slp_inline' performs vortex breeding by nudging sea level pressure toward
!! single point observations
!>@details This is a tropical cyclone 'bogusing' routine that currently only works
!! for hydrostatic dynamics.
!>@note Note: conserve water mass, geopotential, and momentum at the expense of dry air mass
 subroutine breed_slp_inline(nstep, dt, npz, ak, bk, phis, pe, pk, peln, pkz, delp, u, v, pt, q, nwat,   &
                             zvir, gridstruct, ks, domain_local, bd, hydrostatic)
! Input
      integer, intent(in):: nstep, npz, nwat, ks
      real, intent(in):: dt       !< (small) time step in seconds
      real, intent(in):: zvir
      real, intent(in), dimension(npz+1):: ak, bk
      logical, intent(in):: hydrostatic
      type(fv_grid_bounds_type), intent(IN) :: bd
      real, intent(in):: phis(isd:ied,jsd:jed)
      type(domain2d), intent(INOUT) :: domain_local
! Input/Output
      real, intent(inout):: u(isd:ied,jsd:jed+1,npz)
      real, intent(inout):: v(isd:ied+1,jsd:jed,npz)
      real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp, pt
      real, intent(inout)::q(isd:ied,jsd:jed,npz,*)

      real, intent(inout):: pk(is:ie,js:je, npz+1)          !< pe**kappa
      real, intent(inout):: pe(is-1:ie+1, npz+1,js-1:je+1)  !< edge pressure (pascal)
      real, intent(inout):: pkz(is:ie,js:je,npz)
      real, intent(out):: peln(is:ie,npz+1,js:je)           !< ln(pe)

      type(fv_grid_type), target :: gridstruct
! local
      type(time_type):: time
      real:: ps(is:ie,js:je)
      real:: dist(is:ie,js:je)
      real::   tm(is:ie,js:je)
      real::  slp(is:ie,js:je)
      real(kind=R_GRID):: pos(2)
      real:: slp_o         ! sea-level pressure (Pa)
      real:: w10_o, p_env
      real:: r_vor
      real:: relx0, relx, f1, pbreed, pbtop, delp0, dp0
      real:: ratio, p_count, p_sum, a_sum, mass_sink, delps
      real:: p_lo, p_hi, tau_vt, mslp0
      real:: split_time, fac, pdep, r2, r3
      integer year, month, day, hour, minute, second
      integer n, i, j, k, iq, k0

      real, pointer :: area(:,:)
      real(kind=R_GRID), pointer :: agrid(:,:,:)
#ifdef MULTI_GASES
      real :: kappax(is:ie,js:je,npz)
#endif

#if defined (BYPASS_BREED_SLP_INLINE)
    peln = 0.0 ! to silence compiler warning. A dummy argument with an explicit INTENT(OUT) declaration is not given an explicit value.
    call mpp_error(fatal, "breed_slp_inline routine has been disabled")
#else
      if ( forecast_mode ) return

      agrid => gridstruct%agrid_64
      area  => gridstruct%area

    if ( nstorms==0 ) then
         if(master) write(*,*) 'NO TC data to process'
         return
    endif

   if ( k_breed==0 ) k_breed = max(1, ks)


! Advance (local) time
    call get_date(fv_time, year, month, day, hour, minute, second)
    if ( year /= year_track_data ) then
        if (master) write(*,*) 'Warning: The year in storm track data is not the same as model year'
        return
     endif
    time = fv_time   ! fv_time is the time at past time step (set in fv_diag)
    split_time = calday(year, month, day, hour, minute, second) + dt*real(nstep)/86400.

    elapsed_time = elapsed_time + dt
    if ( elapsed_time > nudged_time + 0.1 ) then
         if(print_end_breed)  then
            print_end_breed = .false.
            if (master) write(*,*) '*** Vortext Breeding Ended at', day, hour, minute, second
         endif
         return        !  time to return to forecast mode
    endif

!$OMP parallel do default(none) shared(is,ie,js,je,npz,ps,ak,delp,tm,pkz,pt)
    do j=js,je
! ---- Compute ps
       do i=is,ie
          ps(i,j) = ak(1)
       enddo
       do k=1,npz
          do i=is,ie
             ps(i,j) = ps(i,j) + delp(i,j,k)
          enddo
       enddo
! Compute lowest layer air temperature:
       do i=is,ie
          tm(i,j) = pkz(i,j,npz)*pt(i,j,npz) / cp_air        ! virtual temp
       enddo
    enddo
!   call prt_maxmin('TM', tm, is, ie, js, je, 0, 1, 1.)

!$OMP parallel do default(none) shared(k_breed,npz,conserve_mom,is,ie,js,je,u,v,delp, &
!$OMP                                  conserve_hgt,hydrostatic,pt,pkz,q,nwat)
    do k=k_breed+1,npz

       if ( conserve_mom ) then
       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) * (delp(i,j-1,k)+delp(i,j,k))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) * (delp(i-1,j,k)+delp(i,j,k))
          enddo
       enddo
       endif

       if ( conserve_hgt .and. hydrostatic ) then
       do j=js,je
          do i=is,ie
! Conserve total enthalpy
             pt(i,j,k) = pt(i,j,k)*pkz(i,j,k)*delp(i,j,k)
          enddo
       enddo
       endif

! Convert tracer moist mixing ratio to mass
       do iq=1,nwat
          do j=js,je
             do i=is,ie
                q(i,j,k,iq) = q(i,j,k,iq) * delp(i,j,k)
             enddo
          enddo
       enddo

    enddo

    do 5000 n=1,nstorms      ! loop through all storms

      if ( nobs_tc(n) < min_nobs ) goto 5000

! Check min MSLP
      mslp0 = 1013.E2
      do i=1,nobs_tc(n)
         mslp0 = min( mslp0, mslp_obs(i,n) )
      enddo
      if ( mslp0 > min_mslp ) goto 5000

!----------------------------------------
! Obtain slp observation
!----------------------------------------
      call get_slp_obs(time, nobs_tc(n), x_obs(1,n), y_obs(1,n), wind_obs(1,n),  mslp_obs(1,n), mslp_out(1,n), rad_out(1,n),   &
                       time_tc(1,n), pos(1), pos(2), w10_o, slp_o, r_vor, p_env, stime=split_time, fact=fac)

      if ( slp_o<87500. .or. slp_o>slp_env .or. abs(pos(2))*rad2deg>45. ) then
           goto 5000         ! next storm
      endif


      if(nudge_debug .and. master)    &
         write(*,*) 'Vortex breeding for TC:', n, ' slp=',slp_o/100.,pos(1)*rad2deg,pos(2)*rad2deg

! Determine pbtop (top pressure of vortex breeding)
   k0 = k_breed

   if ( use_high_top ) then
#ifdef HIGH_TEST
       if ( slp_o > 1000.E2 ) then
            pbtop = 850.E2
       else
!           pbtop = max(200.E2, 850.E2-20.*(1000.E2-slp_o))
! mp87:
            pbtop = max(100.E2, 850.E2-25.*(1000.E2-slp_o))
#else
       if ( slp_o > 1000.E2 ) then
            pbtop = 750.E2
       else
            pbtop = max(100.E2, 750.E2-20.*(1000.E2-slp_o))
#endif
       endif
   else
! Lower top for vrotex breeding
      if ( slp_o > 1000.E2 ) then
           pbtop = 900.E2
      else
           pbtop = max(500.E2, 900.E2-5.*(1000.E2-slp_o))  ! mp48
      endif
   endif

    do k=1,npz
       pbreed = ak(k) + bk(k)*1.E5
       if ( pbreed>pbtop ) then
            k0 = k
            exit
       endif
    enddo
    k0 = max(k0, k_breed)

      do j=js, je
         do i=is, ie
            dist(i,j) = great_circle_dist( pos, agrid(i,j,1:2), radius)
         enddo
      enddo

      call compute_slp(is, ie, js, je, tm, ps, phis(is:ie,js:je), slp)

    if ( r_vor < 30.E3 .or. p_env<900.E2 ) then

! Compute r_vor & p_env
         r_vor = r_min + (slp_env-slp_o)/25.E2*r_inc

123   continue
      p_count = 0.
        p_sum = 0.
        a_sum = 0.
      do j=js, je
         do i=is, ie
            if( dist(i,j)<(r_vor+del_r) .and. dist(i,j)>r_vor .and. phis(i,j)<250.*grav ) then
                p_count = p_count + 1.
                  p_sum = p_sum + slp(i,j)*area(i,j)
                  a_sum = a_sum + area(i,j)
            endif
         enddo
      enddo

      call mp_reduce_sum(p_count)

      if ( p_count<32. ) then
           if(nudge_debug .and. master) write(*,*) p_count, 'Skipping obs: too few p_count'
           goto 5000
      endif

      call mp_reduce_sum(p_sum)
      call mp_reduce_sum(a_sum)
      p_env = p_sum / a_sum

      if(nudge_debug .and. master) write(*,*) 'Environmental SLP=', p_env/100., ' computed radius=', r_vor/1.E3

      if ( p_env>1015.E2 .or. p_env<920.E2 ) then
         if( nudge_debug ) then
            if(master)  write(*,*) 'Environmental SLP out of bound; skipping obs. p_count=', p_count, p_sum
            call prt_maxmin('SLP_breeding', slp, is, ie, js, je, 0, 1, 0.01)
         endif
         goto 5000
      endif

    endif


      if ( p_env < max(pre0_env, slp_o + 200.0) ) then

         r_vor = r_vor + 25.E3

         if(nudge_debug .and. master) then
            write(*,*) 'Computed environmental SLP too low'
            write(*,*) ' ', p_env/100., slp_o/100.,pos(1)*rad2deg, pos(2)*rad2deg
         endif

         if ( slp_o > 1003.E2 .and.  r_vor > 1000.E3 ) then
!             if(master) write(*,*) 'Failed to size the Vortex for the weak storm'
              goto 5000
         endif

         if ( r_vor < 1250.E3 ) goto 123

!        if(master) write(*,*) 'Failed to size the Vortex; skipping this storm'
         goto 5000

      endif

      tau_vt = tau_vt_slp * (1. + (960.E2-slp_o)/100.E2 )
      tau_vt = max(abs(dt), tau_vt)

      if ( do_adiabatic_init ) then
           relx0  = min(1., 2.*abs(dt)/tau_vt)
      else
           relx0  = min(1., abs(dt)/tau_vt)
      endif

      mass_sink = 0.

!$OMP parallel do default(none) shared(is,ie,js,je,dist,r_vor,phis,p_env,slp_o,r_hi,r_lo, &
!$OMP                                  ps,tm,relx0,tau_vt_rad,dps_min,slp,ak,k0,delp,npz,area) &
!$OMP                           private(f1, p_hi, p_lo, relx, delps, pbreed, mass_sink, dp0, pdep)
      do j=js, je
         do i=is, ie
            if( dist(i,j) < r_vor .and. phis(i,j)<250.*grav ) then
                f1 = dist(i,j)/r_vor
! Compute p_obs: assuming local radial distributions of slp are Gaussian
                p_hi = p_env - (p_env-slp_o) * exp( -r_hi*f1**2 )    ! upper bound
                p_lo = p_env - (p_env-slp_o) * exp( -r_lo*f1**2 )    ! lower bound

                if ( ps(i,j) > p_hi .and. tm(i,j) < tm_max ) then
!                                         do nothing if lowest layer is too hot
! Under-development:
                      relx = relx0*exp( -tau_vt_rad*f1**2 )
                     delps = relx*(ps(i,j) - p_hi) * min( (tm_max-tm(i,j))/10., 1.)
! After mp115
!                    delps = relx*(ps(i,j) - p_hi) * min( (tm_max-tm(i,j))/5., 1.)
                                                    ! Cap the increment to prevent overheating
                                                    ! Note: ps is used here to prevent
                                                    !       over deepening over terrain
                     delps = min(delps, dps_min)
                elseif ( slp(i,j) < p_lo ) then
! Over-development:
                      relx = max(0.5, relx0)
                     delps = relx*(slp(i,j) - p_lo)  ! Note: slp is used here
                else
                     goto 400        ! do nothing; proceed to next storm
                endif

#ifdef SIM_TEST
                pbreed = ak(1)
                do k=1,k0
                   pbreed = pbreed + delp(i,j,k)
                enddo
                f1 = 1. - delps/(ps(i,j)-pbreed)
                do k=k0+1,npz
                   delp(i,j,k) = delp(i,j,k)*f1
                enddo
                mass_sink = mass_sink + delps*area(i,j)
#else
                if ( delps > 0. ) then
                      pbreed = ak(1)
                      do k=1,k0
                         pbreed = pbreed + delp(i,j,k)
                      enddo
                      f1 = 1. - delps/(ps(i,j)-pbreed)
                      do k=k0+1,npz
                         delp(i,j,k) = delp(i,j,k)*f1
                      enddo
                      mass_sink = mass_sink + delps*area(i,j)
                else
                      dp0 = abs(delps)
                      do k=npz,k0+1,-1
                         if ( abs(delps) < 1. ) then
                              delp(i,j,k) = delp(i,j,k) - delps
                              mass_sink = mass_sink + delps*area(i,j)
                              go to 400
                         else
                              pdep = max(1.0, min(abs(0.4*delps), 0.2*dp0,  0.02*delp(i,j,k)))
                              pdep = - min(pdep, abs(delps))
                              delp(i,j,k) = delp(i,j,k) - pdep
                              delps = delps - pdep
                              mass_sink = mass_sink + pdep*area(i,j)
                         endif
                      enddo
                endif
#endif

            endif
400     continue
        enddo        ! end i-loop
      enddo        ! end j-loop

      call mp_reduce_sum(mass_sink)
      if ( abs(mass_sink)<1.E-40 ) goto 5000

      r2 = r_vor + del_r
      r3 = min( 4.*r_vor, max(2.*r_vor, 2500.E3) ) + del_r

      p_sum = 0.
      do j=js, je
         do i=is, ie
            if( dist(i,j)<r3 .and. dist(i,j)>r2 ) then
                p_sum = p_sum + area(i,j)
            endif
         enddo
      enddo

      call mp_reduce_sum(p_sum)
      mass_sink = mass_sink / p_sum ! mean delta pressure to be added back to the environment to conserve mass
      if(master .and. nudge_debug) write(*,*) 'TC#',n, 'Mass tele-ported (Pa)=', mass_sink

!$OMP parallel do default(none) shared(is,ie,js,je,dist,r3,r2,ak,k_breed,delp,ps,mass_sink,npz) &
!$OMP             private(pbreed, f1)
      do j=js, je
         do i=is, ie
            if( dist(i,j)<r3 .and. dist(i,j)>r2 ) then
                pbreed = ak(1)
                do k=1,k_breed
                   pbreed = pbreed + delp(i,j,k)
                enddo
                f1 = 1. + mass_sink/(ps(i,j)-pbreed)
                do k=k_breed+1,npz
                   delp(i,j,k) = delp(i,j,k)*f1
                enddo
            endif
         enddo
      enddo

! ---- re-compute ps
      do j=js,je
         do i=is,ie
            ps(i,j) = ak(1)
         enddo
         do k=1,npz
            do i=is,ie
               ps(i,j) = ps(i,j) + delp(i,j,k)
            enddo
         enddo
      enddo

5000 continue

!--------------------------
! Update delp halo regions:
!--------------------------
    call mpp_update_domains(delp, domain_local, complete=.true.)

!$OMP parallel do default(none) shared(is,ie,js,je,npz,pe,ak,delp,k_breed,peln,pk)
    do j=js-1,je+1
       do i=is-1,ie+1
          pe(i,1,j) = ak(1)
       enddo
       do k=2,npz+1
          do i=is-1,ie+1
             pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
          enddo
       enddo
    enddo

#ifdef MULTI_GASES
!$OMP parallel do default(none) shared(is,ie,js,je,npz,k_breed,kappax,num_gas)
    do k=k_breed,npz
      do j=js,je
         do i=is,ie
            kappax(i,j,k)= virqd(q(i,j,k,1:num_gas))/vicpqd(q(i,j,k,1:num_gas))
         enddo
      enddo
    enddo
#endif
    do k=k_breed+1,npz+1
      do j=js,je
         do i=is,ie
            peln(i,k,j) = log(pe(i,k,j))
              pk(i,j,k) = pe(i,k,j)**kappa
         enddo
      enddo
    enddo

!$OMP parallel do default(none) shared(k_breed,npz,is,ie,js,je,conserve_mom,u,v,delp, &
#ifdef MULTI_GASES
!$OMP                                  kappax,                                        &
#endif
!$OMP                                  conserve_hgt,hydrostatic,pkz,pk,pt,peln)
    do k=k_breed+1,npz

       if ( conserve_mom ) then
       do j=js,je+1
          do i=is,ie
             u(i,j,k) = u(i,j,k) / (delp(i,j-1,k)+delp(i,j,k))
          enddo
       enddo
       do j=js,je
          do i=is,ie+1
             v(i,j,k) = v(i,j,k) / (delp(i-1,j,k)+delp(i,j,k))
          enddo
       enddo
       endif

       if ( conserve_hgt .and. hydrostatic ) then
       do j=js,je
          do i=is,ie
! Conserve total enthalpy (static energy)
            pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
#ifdef MULTI_GASES
            pkz(i,j,k) = exp(kappax(i,j,k)*log(pkz(i,j,k)))
#endif
             pt(i,j,k) = pt(i,j,k) / (pkz(i,j,k)*delp(i,j,k))
          enddo
       enddo
! pkz last used
       endif
    enddo

! Convert tracer mass back to moist mixing ratio

!$OMP parallel do default(none) shared(nwat,k_breed,npz,is,ie,js,je,q,delp)
    do iq=1,nwat
       do k=k_breed+1,npz
          do j=js,je
             do i=is,ie
                q(i,j,k,iq) = q(i,j,k,iq) / delp(i,j,k)
             enddo
          enddo
       enddo
    enddo

    if(hydrostatic) call mpp_update_domains(pt, domain_local, complete=.true.)

    nullify(agrid)
    nullify(area)
#endif

  end subroutine breed_slp_inline

!>@brief The subroutine 'breed_srf_w10' performs vortex breeding by nudging 10-m winds.
!>@details This is the inline version.
 subroutine breed_srf_w10(time, dt, npz, ak, bk, ps, phis, slp, delp, u, v, gridstruct)
      type(time_type), intent(in):: time
      integer, intent(in):: npz
      real, intent(in):: dt       !< time step in seconds
      real, intent(in), dimension(npz+1):: ak, bk
      real, intent(in):: phis(isd:ied,jsd:jed)
      real, intent(in)::   ps(isd:ied,jsd:jed)
      real, intent(in), dimension(is:ie,js:je):: slp
      type(fv_grid_type), intent(IN), target :: gridstruct
      real, intent(inout):: u(isd:ied,jsd:jed+1,npz)
      real, intent(inout):: v(isd:ied+1,jsd:jed,npz)
! Input/Output
      real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp
! local
      real, dimension(is:ie,js:je):: us, vs
      real wu(is:ie,  js:je+1)
      real wv(is:ie+1,js:je)
      real u1(is:ie), v1(is:ie)
      real:: dist(is:ie,js:je), wind(is:ie,js:je)
      real(kind=R_GRID):: pos(2)
      real:: slp_o         !< sea-level pressure (Pa)
      real:: w10_o, p_env
      real:: r_vor, pc, p_count
      real:: r_max, speed, ut, vt, speed_local        ! tangent wind speed
      real:: u_bg, v_bg, mass, t_mass
      real:: relx0, relx, f1
      real:: z0, mslp0
      real:: zz = 35.           !< mid-layer height at the lowest model level
      real:: wind_fac           !< Computed ( ~ 1.2)
      integer n, i, j

      ! Pointers
      real, pointer, dimension(:,:) :: dx, dy, rdxa, rdya, a11, a21, a12, a22, area
      real(kind=R_GRID), pointer, dimension(:,:,:) :: agrid, vlon, vlat

      dx    => gridstruct%dx
      dy    => gridstruct%dy
      rdxa  => gridstruct%rdxa
      rdya  => gridstruct%rdya
      a11   => gridstruct%a11
      a21   => gridstruct%a21
      a12   => gridstruct%a12
      a22   => gridstruct%a22
      area  => gridstruct%area
      agrid => gridstruct%agrid_64
      vlon  => gridstruct%vlon
      vlat  => gridstruct%vlat


    if ( nstorms==0 ) then
         if(master) write(*,*) 'NO TC data to process'
         return
    endif

! Compute lat-lon winds on A grid
!$OMP parallel do default(none) shared(is,ie,js,je,wu,u,npz,dx)
    do j=js,je+1
       do i=is,ie
          wu(i,j) = u(i,j,npz)*dx(i,j)
       enddo
    enddo
!$OMP parallel do default(none) shared(is,ie,js,je,wv,v,npz,dy)
    do j=js,je
       do i=is,ie+1
          wv(i,j) = v(i,j,npz)*dy(i,j)
       enddo
    enddo

!$OMP parallel do default(none) shared(is,ie,js,je,wu,rdxa,wv,rdya,us,vs,a11,a12,a21,a22,wind) &
!$OMP                          private(u1,v1)
    do j=js, je
       do i=is, ie
! Co-variant to Co-variant "vorticity-conserving" interpolation
             u1(i) = (wu(i,j) + wu(i,j+1)) * rdxa(i,j)
             v1(i) = (wv(i,j) + wv(i+1,j)) * rdya(i,j)
! Cubed (cell center co-variant winds) to lat-lon:
             us(i,j) = a11(i,j)*u1(i) + a12(i,j)*v1(i)
             vs(i,j) = a21(i,j)*u1(i) + a22(i,j)*v1(i)
! Wind speed
          wind(i,j) = sqrt( us(i,j)**2 + vs(i,j)**2 )
       enddo
    enddo

    relx0  = min(1., dt/tau_vt_wind)

    do 3000 n=1,nstorms  ! loop through all storms

      if ( nobs_tc(n) < min_nobs ) goto 3000

! Check min MSLP
      mslp0 = 1013.E2
      do i=1,nobs_tc(n)
         mslp0 = min( mslp0, mslp_obs(i,n) )
      enddo
      if ( mslp0 > min_mslp ) goto 3000

!----------------------------------------
! Obtain slp observation
!----------------------------------------
      call get_slp_obs(time, nobs_tc(n), x_obs(1,n), y_obs(1,n), wind_obs(1,n),  mslp_obs(1,n), mslp_out(1,n), rad_out(1,n),   &
                       time_tc(1,n), pos(1), pos(2), w10_o, slp_o, r_vor, p_env)

      if ( slp_o<90000. .or. slp_o>slp_env .or. abs(pos(2))*rad2deg>35. ) goto 3000         ! next storm


      do j=js, je
         do i=is, ie
            dist(i,j) = great_circle_dist( pos, agrid(i,j,1:2), radius )
         enddo
      enddo

      r_vor = r_min + (slp_env-slp_o)/25.E2*r_inc

!----------------------------------------------------
! * Find model's SLP center nearest to the observation
! * Find maximum wind speed at the lowest model level

     speed_local = 0.
           r_max = -999.
              pc = 1013.E2
     p_count = 0.
     do j=js, je
        do i=is, ie
           if( dist(i,j) < r_vor ) then

! Counting the "land" points:
               if ( dist(i,j)<0.5*r_vor .and. abs(phis(i,j))>2.*grav )  p_count = p_count + 1.

               pc = min(pc, slp(i,j))

               if ( speed_local < wind(i,j) ) then
                    speed_local = wind(i,j)
                    r_max = dist(i,j)
               endif

           endif
        enddo
     enddo

     call mp_reduce_sum(p_count)
     if ( p_count>32 ) goto 3000  ! over/near rough land

     if ( w10_o < 0. ) then   ! 10-m wind obs is not available
! Uses Atkinson_Holliday wind-pressure correlation
          w10_o = 3.446778 * (1010.-slp_o/100.)**0.644
     endif

     speed = speed_local
     call mp_reduce_max(speed)     ! global max wind (near storm)
     call mp_reduce_min(pc)

    if ( speed_local < speed ) then
         r_max = -999.
    endif
    call mp_reduce_max(r_max)
    if( r_max<0. ) call mpp_error(FATAL,'==> Error in r_max')

! ---------------------------------------------------
! Determine surface wind speed and radius for nudging
! ---------------------------------------------------

! Compute surface roughness z0 from w10, based on Eq (4) & (5) from Moon et al. 2007
     if ( w10_o > 12.5 ) then
          z0 = (0.085*w10_o - 0.58) * 1.E-3
! z0 (w10=40) = 2.82E-3
!         z0 = min( z0, 2.82E-3 )

     else
          z0 = 0.0185/grav*(0.001*w10_o**2 + 0.028*w10_o)**2
     endif

! lowest layer height: zz

     wind_fac = log(zz/z0) / log(10./z0)
     if( nudge_debug .and. master ) write(*,*) 'Wind adjustment factor=', wind_fac
     if( wind_fac<1. ) call mpp_error(FATAL,'==> Error in wind_fac')

     if ( pc < slp_o ) then
!--
!         The storm in the model is over developed
!         if ( (pc+3.0E2)>slp_o .or. speed <= w10_o ) go to 3000    ! next storm
!--
! using radius (r_max) as dtermined above;
! What if model's pressure center is very far from the observed?
! using obs wind
          speed = wind_fac*w10_o
     else
!         The storm in the model is under developed; using max wind
          speed = max(wind_fac*w10_o, speed)
          if ( pc>1009.E2 )  r_max = 0.5 * r_vor
     endif

! More adjustment

! Some bounds on the radius of maximum wind:
     r_max = max(2.5*grid_size, r_max)      ! at least 2.5X the grid size
     r_max = min(0.75*r_vor, r_max)

     t_mass = 0.
     u_bg = 0.
     v_bg = 0.

     if ( add_bg_wind ) then
       p_count = 0.
       do j=js, je
          do i=is, ie
           if( dist(i,j) <= min(r_vor,2.*r_fac*r_max) .and. phis(i,j)<2.*grav ) then
               mass = area(i,j)*delp(i,j,npz)
!-- using model winds ----------------------------------
               u_bg = u_bg + us(i,j)*mass
               v_bg = v_bg + vs(i,j)*mass
!-------------------------------------------------------
               t_mass = t_mass + mass
               p_count = p_count + 1.
           endif
          enddo
       enddo
       call mp_reduce_sum(p_count)
       if ( p_count<16. ) go to 3000

       call mp_reduce_sum(t_mass)
       call mp_reduce_sum(u_bg)
       call mp_reduce_sum(v_bg)
       u_bg = u_bg / t_mass
       v_bg = v_bg / t_mass
       if ( nudge_debug .and. master ) write(*,*) pos(2)*rad2deg, 'vortex bg wind=', u_bg, v_bg
     endif

      relx = relx0
! Nudge wind in the "inner core":
#ifdef TTT
      do j=js, je
         do i=is, ie
            if( dist(i,j) <= min(r_vor, r_fac*r_max) .and. phis(i,j)<2.*grav ) then
                f1 = dist(i,j)/r_max
                relx = relx0*exp( -tau_vt_rad*f1**2 )
                if( dist(i,j)<=r_max ) then
                    speed_local = speed * f1
                else
                    speed_local = speed / f1**0.75
                endif
                call tangent_wind(vlon(i,j,1:3), vlat(i,j,1:3), speed_local, pos, agrid(i,j,1:2), ut, vt)
                ut = ut + u_bg
                vt = vt + v_bg
! Update:
                us(i,j) = relx*(ut-us(i,j))
                vs(i,j) = relx*(vt-vs(i,j))
            endif
400     continue
        enddo        ! end i-loop
      enddo        ! end j-loop
#else
      do j=js, je
         do i=is, ie
            if( dist(i,j) <= min(r_vor, r_fac*r_max) .and. phis(i,j)<2.*grav ) then
                f1 = dist(i,j)/r_max
                relx = relx0*exp( -tau_vt_rad*f1**2 )
                if( dist(i,j)<=r_max ) then
                    speed_local = speed * f1
                else
                    speed_local = speed / f1**0.75
                endif
                call tangent_wind(vlon(i,j,1:3), vlat(i,j,1:3), speed_local, pos, agrid(i,j,1:2), ut, vt)
                ut = ut + u_bg
                vt = vt + v_bg
! Update:
                us(i,j) = relx*(ut-us(i,j))
                vs(i,j) = relx*(vt-vs(i,j))
            endif
400     continue
        enddo        ! end i-loop
      enddo        ! end j-loop
#endif

3000 continue

  end subroutine breed_srf_w10

!>@brief The subroutine 'breed_srf_winds' performs vortex breeding by nudging 1-m winds.
 subroutine breed_srf_winds(time, dt, npz, u_obs, v_obs, ak, bk, ps, phis, delp, ua, va, u_dt, v_dt, pt, q, nwat, zvir, gridstruct)
! Input
      type(time_type), intent(in):: time
      integer, intent(in):: npz, nwat
      real, intent(in):: dt       !< time step in seconds
      real, intent(in):: zvir
      real, intent(in), dimension(npz+1):: ak, bk
      real, intent(in), dimension(isd:ied,jsd:jed):: phis, ps
      real, intent(in), dimension(is:ie,js:je,npz):: u_obs, v_obs
      type(fv_grid_type), intent(IN), target :: gridstruct
! Input/Output
      real, intent(inout), dimension(isd:ied,jsd:jed,npz):: delp, pt, ua, va, u_dt, v_dt
      real, intent(inout)::q(isd:ied,jsd:jed,npz,nwat)
! local
      real:: dist(is:ie,js:je), wind(is:ie,js:je)
      real::  slp(is:ie,js:je)
      real(kind=R_GRID):: pos(2)
      real:: slp_o         !< sea-level pressure (Pa)
      real:: w10_o, p_env
      real:: r_vor, pc, p_count
      real:: r_max, speed, ut, vt, speed_local        ! tangent wind speed
      real:: u_bg, v_bg, mass, t_mass
      real:: relx0, relx, f1, rdt
      real:: z0, mslp0
      real:: zz = 35.           !< mid-layer height at the lowest model level
      real:: wind_fac           !< Computed ( ~ 1.2)
      integer n, i, j, k, iq

      real, pointer :: area(:,:)
      real(kind=R_GRID), pointer :: vlon(:,:,:), vlat(:,:,:), agrid(:,:,:)

      area => gridstruct%area
      vlon => gridstruct%vlon
      vlat => gridstruct%vlat
      agrid => gridstruct%agrid_64

    if ( nstorms==0 ) then
         if(master) write(*,*) 'NO TC data to process'
         return
    endif

       rdt = 1./dt
    relx0  = min(1., dt/tau_vt_wind)

!$OMP parallel do default(none) shared(is,ie,js,je,slp,ps,phis,pt,npz,wind,ua,va)
    do j=js, je
       do i=is, ie
           slp(i,j) = ps(i,j)*exp(phis(i,j)/(rdgas*(pt(i,j,npz)+3.25E-3*phis(i,j)/grav)))
          wind(i,j) = sqrt( ua(i,j,npz)**2 + va(i,j,npz)**2 )
       enddo
    enddo

!!!!!$OMP parallel do default(none) private(pos, w10_o, slp_o, r_vor, p_env)
    do 3000 n=1,nstorms  ! loop through all storms

      if ( nobs_tc(n) < min_nobs ) goto 3000

! Check min MSLP
      mslp0 = 1013.E2
      do i=1,nobs_tc(n)
         mslp0 = min( mslp0, mslp_obs(i,n) )
      enddo
      if ( mslp0 > min_mslp ) goto 3000

!----------------------------------------
! Obtain slp observation
!----------------------------------------
      call get_slp_obs(time, nobs_tc(n), x_obs(1,n), y_obs(1,n), wind_obs(1,n),  mslp_obs(1,n), mslp_out(1,n), rad_out(1,n),   &
                       time_tc(1,n), pos(1), pos(2), w10_o, slp_o, r_vor, p_env)

      if ( slp_o<90000. .or. slp_o>slp_env .or. abs(pos(2))*rad2deg>35. ) goto 3000         ! next storm


      do j=js, je
         do i=is, ie
            dist(i,j) = great_circle_dist( pos, agrid(i,j,1:2), radius )
         enddo
      enddo

      r_vor = r_min + (slp_env-slp_o)/25.E2*r_inc

!----------------------------------------------------

! * Find model's SLP center nearest to the observation
! * Find maximum wind speed at the lowest model level

     speed_local = 0.
           r_max = -999.
              pc = 1013.E2
     p_count = 0.
     do j=js, je
        do i=is, ie
           if( dist(i,j) < r_vor ) then

! Counting the "land" points:
               if ( dist(i,j)<0.5*r_vor .and. abs(phis(i,j))>2.*grav )  p_count = p_count + 1.

               pc = min(pc, slp(i,j))

               if ( speed_local < wind(i,j) ) then
                    speed_local = wind(i,j)
                    r_max = dist(i,j)
               endif

           endif
        enddo
     enddo

     call mp_reduce_sum(p_count)
     if ( p_count>32 ) goto 3000  ! over/near rough land

     if ( w10_o < 0. ) then   ! 10-m wind obs is not available
! Uses Atkinson_Holliday wind-pressure correlation
          w10_o = 3.446778 * (1010.-slp_o/100.)**0.644
     endif

     speed = speed_local
     call mp_reduce_max(speed)     ! global max wind (near storm)
     call mp_reduce_min(pc)

    if ( speed_local < speed ) then
         r_max = -999.
    endif
    call mp_reduce_max(r_max)
    if( r_max<0. ) call mpp_error(FATAL,'==> Error in r_max')

! ---------------------------------------------------
! Determine surface wind speed and radius for nudging
! ---------------------------------------------------

! Compute surface roughness z0 from w10, based on Eq (4) & (5) from Moon et al. 2007
     if ( w10_o > 12.5 ) then
          z0 = (0.085*w10_o - 0.58) * 1.E-3
! z0 (w10=40) = 2.82E-3
!         z0 = min( z0, 2.82E-3 )

     else
          z0 = 0.0185/grav*(0.001*w10_o**2 + 0.028*w10_o)**2
     endif

! lowest layer height: zz

     wind_fac = log(zz/z0) / log(10./z0)
     if( nudge_debug .and. master ) write(*,*) 'Wind adjustment factor=', wind_fac
     if( wind_fac<1. ) call mpp_error(FATAL,'==> Error in wind_fac')

     if ( pc < slp_o ) then
!--
!         The storm in the model is over developed
!         if ( (pc+3.0E2)>slp_o .or. speed <= w10_o ) go to 3000    ! next storm
!--
! using radius (r_max) as dtermined above;
! What if model's pressure center is very far from the observed?
! using obs wind
          speed = wind_fac*w10_o
     else
!         The storm in the model is under developed; using max wind
          speed = max(wind_fac*w10_o, speed)
          if ( pc>1009.E2 )  r_max = 0.5 * r_vor
     endif

! More adjustment

! Some bounds on the radius of maximum wind:
     r_max = max(2.5*grid_size, r_max)      ! at least 2.5X the grid size
     r_max = min(0.75*r_vor, r_max)

     t_mass = 0.
     u_bg = 0.
     v_bg = 0.

     if ( add_bg_wind ) then
       p_count = 0.
       do j=js, je
          do i=is, ie
           if( dist(i,j) <= min(r_vor,2.*r_fac*r_max) .and. phis(i,j)<2.*grav ) then
               mass = area(i,j)*delp(i,j,npz)
!-- using model winds ----------------------------------
!              u_bg = u_bg + ua(i,j,npz)*mass
!              v_bg = v_bg + va(i,j,npz)*mass
!-------------------------------------------------------
! Using analysis winds
               u_bg = u_bg + u_obs(i,j,npz)*mass
               v_bg = v_bg + v_obs(i,j,npz)*mass
               t_mass = t_mass + mass
               p_count = p_count + 1.
           endif
          enddo
       enddo
       call mp_reduce_sum(p_count)
       if ( p_count<16. ) go to 3000

       call mp_reduce_sum(t_mass)
       call mp_reduce_sum(u_bg)
       call mp_reduce_sum(v_bg)
       u_bg = u_bg / t_mass
       v_bg = v_bg / t_mass
!      if ( master ) write(*,*) pos(2)*rad2deg, 'vortex bg wind=', u_bg, v_bg
     endif

     relx = relx0
     k = npz                 ! lowest layer only
! Nudge wind in the "inner core":
      do j=js, je
         do i=is, ie
            if( dist(i,j) <= min(r_vor, r_fac*r_max) .and. phis(i,j)<2.*grav ) then
                f1 = dist(i,j)/r_max
                relx = relx0*exp( -tau_vt_rad*f1**2 )
                if( dist(i,j)<=r_max ) then
                    speed_local = speed * f1
                else
                    speed_local = speed / f1**0.75
                endif
                call tangent_wind(vlon(i,j,1:3), vlat(i,j,1:3), speed_local, pos, agrid(i,j,1:2), ut, vt)
                ut = ut + u_bg
                vt = vt + v_bg
                u_dt(i,j,k) = u_dt(i,j,k) + relx*(ut-ua(i,j,k)) * rdt
                v_dt(i,j,k) = v_dt(i,j,k) + relx*(vt-va(i,j,k)) * rdt
! Update:
                ua(i,j,k) = ua(i,j,k) + relx*(ut-ua(i,j,k))
                va(i,j,k) = va(i,j,k) + relx*(vt-va(i,j,k))
            endif
400     continue
        enddo        ! end i-loop
      enddo        ! end j-loop

3000 continue

  end subroutine breed_srf_winds

  subroutine tangent_wind ( elon, elat, speed, po, pp, ut, vt )
  real, intent(in):: speed
  real(kind=R_GRID), intent(in):: po(2), pp(2)
  real(kind=R_GRID), intent(in):: elon(3), elat(3)
  real, intent(out):: ut, vt
! local
  real(kind=R_GRID):: e1(3), eo(3), ep(3), op(3)

  call latlon2xyz(po, eo)
  call latlon2xyz(pp, ep)

  op(:) = ep(:) - eo(:)
  call normalize_vect( op )

  call vect_cross(e1, ep, eo)

  ut = speed * (e1(1)*elon(1) + e1(2)*elon(2) + e1(3)*elon(3))
  vt = speed * (e1(1)*elat(1) + e1(2)*elat(2) + e1(3)*elat(3))

! SH:
  if ( po(2) < 0. ) then
       ut = -ut
       vt = -vt
  endif

  end subroutine tangent_wind

  subroutine get_slp_obs(time, nobs, lon_obs, lat_obs, w10, mslp, slp_out, r_out, time_obs,    &
                         x_o, y_o, w10_o, slp_o, r_vor, p_vor, stime, fact)
! Input
    type(time_type), intent(in):: time
    integer, intent(in)::  nobs   !< number of observations in this particular storm
    real(KIND=4), intent(in)::  lon_obs(nobs)
    real(KIND=4), intent(in)::  lat_obs(nobs)
    real(KIND=4), intent(in)::      w10(nobs)        !< observed 10-m wind speed
    real(KIND=4), intent(in)::     mslp(nobs)        !< observed SLP in Pa
    real(KIND=4), intent(in)::  slp_out(nobs)        !< slp at r_out
    real(KIND=4), intent(in)::    r_out(nobs)
    real(KIND=4), intent(in):: time_obs(nobs)
    real, optional, intent(in):: stime
    real, optional, intent(out):: fact
! Output
    real(kind=R_GRID), intent(out):: x_o , y_o      !< position of the storm center
    real, intent(out):: w10_o          !< 10-m wind speed
    real, intent(out):: slp_o          !< Observed sea-level-pressure (Pa)
    real, intent(out):: r_vor, p_vor
! Internal:
    real:: t_thresh
      real(kind=R_GRID):: p1(2), p2(2)
      real time_model
      real(kind=R_GRID) fac
      integer year, month, day, hour, minute, second, n

      t_thresh = 600./86400.  ! unit: days

       w10_o = -100000.
       slp_o = -100000.
         x_o = -100.*pi
         y_o = -100.*pi
       p_vor = -1.E10
       r_vor = -1.E10

   if ( present(stime) ) then
      time_model = stime
   else
      call get_date(time, year, month, day, hour, minute, second)

      if ( year /= year_track_data ) then
           if (master) write(*,*) 'Warning: The year in storm track data is not the same as model year'
           return
      endif

      time_model = calday(year, month, day, hour, minute, second)
!     if(nudge_debug .and. master) write(*,*) 'Model:', time_model, year, month, day, hour, minute, second
   endif

!-------------------------------------------------------------------------------------------
!     if ( time_model <= time_obs(1)  .or.  time_model >= time_obs(nobs) ) then
!          return
!-------------------------------------------------------------------------------------------

      if ( time_model <= (time_obs(1)-t_thresh)  .or.  time_model >= time_obs(nobs) ) return

      if ( time_model <=  time_obs(1) ) then
!--
! This is an attempt to perform vortex breeding several minutes before the first available observation
!--
                 w10_o =     w10(1)
                 slp_o =    mslp(1)
                   x_o = lon_obs(1)
                   y_o = lat_obs(1)
                 if ( present(fact) )  fact = 1.25
      else
           do n=1,nobs-1
             if( time_model >= time_obs(n) .and. time_model <= time_obs(n+1) ) then
                   fac = (time_model-time_obs(n)) / (time_obs(n+1)-time_obs(n))
                 w10_o =     w10(n) + (    w10(n+1)-    w10(n)) * fac
                 slp_o =    mslp(n) + (   mslp(n+1)-   mslp(n)) * fac
! Trajectory interpolation:
! Linear in (lon,lat) space
!                  x_o = lon_obs(n) + (lon_obs(n+1)-lon_obs(n)) * fac
!                  y_o = lat_obs(n) + (lat_obs(n+1)-lat_obs(n)) * fac
                 p1(1) = lon_obs(n);     p1(2) = lat_obs(n)
                 p2(1) = lon_obs(n+1);   p2(2) = lat_obs(n+1)
                 call intp_great_circle(fac, p1, p2, x_o, y_o)
!----------------------------------------------------------------------
                  if ( present(fact) )   fact = 1. + 0.25*cos(fac*2.*pi)
! Additional data from the extended best track
!                if ( slp_out(n)>0. .and. slp_out(n+1)>0. .and. r_out(n)>0. .and. r_out(n+1)>0. ) then
!                     p_vor = slp_out(n) + ( slp_out(n+1) - slp_out(n)) * fac
!                     r_vor =   r_out(n) + (   r_out(n+1) -   r_out(n)) * fac
!                endif
                 return
             endif
           enddo
      endif

  end subroutine get_slp_obs


  subroutine slp_obs_init
  integer:: unit, n, nobs
  character(len=3):: GMT
  character(len=9):: ts_name
  character(len=19):: comment
  integer:: mmddhh, yr, year, month, day, hour, MPH, islp
  integer:: it, i1, i2, p_ring, d_ring
  real:: lon_deg, lat_deg, cald, slp, mps

  nobs_tc(:) = 0
  time_tc(:,:) = 0.
  wind_obs(:,:) = -100000.
  mslp_obs(:,:) = -100000.
  x_obs(:,:) = - 100.*pi
  y_obs(:,:) = - 100.*pi

  mslp_out(:,:) = -1.E10
   rad_out(:,:) = -1.E10

  if( track_file_name == "No_File_specified" ) then
      if(master) write(*,*) 'No TC track file specified'
      return
  else
      unit = get_unit()
      open( unit, file=track_file_name)
  endif

  read(unit, *) year
  if(master) write(*,*) 'Reading TC track data for YEAR=', year

  year_track_data = year

  nstorms = 0
     nobs = 0
    month = 99

  if ( ibtrack ) then

!---------------------------------------------------------------
! The data format is from Ming Zhoa's processed ibTrack datasets
!---------------------------------------------------------------

    read(unit, *) ts_name, nobs, yr, month, day, hour

    if ( yr /= year ) then
         if(master) write(*, *) 'Year inconsistency found !!!'
         call mpp_error(FATAL,'==> Error in reading best track data')
    endif

    do while ( ts_name=='start' )

               nstorms  = nstorms + 1
       nobs_tc(nstorms) = nobs       ! observation count for this storm
       if(nudge_debug .and. master) write(*, *) 'Read Data for TC#', nstorms, nobs

       do it=1, nobs
          read(unit, *) lon_deg, lat_deg, mps, slp, yr, month, day, hour
!         if ( yr /= year ) then
!             if(master) write(*, *) 'Extended to year + 1', yr
!         endif
          cald = calday(yr, month, day, hour, 0, 0)
          time_tc(it,nstorms) = cald
          if(nudge_debug .and. master) write(*, 100) cald, month, day, hour, lon_deg, lat_deg, mps, slp

          wind_obs(it,nstorms) = mps       ! m/s
          mslp_obs(it,nstorms) = 100.*slp
             y_obs(it,nstorms) = lat_deg * deg2rad
             x_obs(it,nstorms) = lon_deg * deg2rad
       enddo

       read(unit, *) ts_name, nobs, yr, month, day, hour
    enddo
100  format(1x, f9.2, 1x, i3, 1x, i2, 1x, i2, 1x, f6.1, 1x, f6.1, 1x, f4.1, 1x, f6.1)

  else

  do while ( month /= 0 )

     read(unit, *) month, day, hour, GMT, lat_deg, lon_deg, MPH, islp, comment

     select case (month)

     case (99)                ! New storm record to start
          nstorms = nstorms + 1
          nobs = 0
          if(master) write(*, *) 'Reading data for TC#', nstorms, comment
     case ( 0)                ! end of record
          if(master) write(*, *) 'End of record reached'
     case default
           nobs = nobs + 1
           cald = calday(year, month, day, hour, 0, 0)
           time_tc(nobs,nstorms) = cald
           nobs_tc(nstorms) = nobs       ! observation count for this storm

          if(master) write(*, 200) nobs, cald,  month, day, hour, GMT, lat_deg, lon_deg, MPH, islp, comment
          mslp_obs(nobs,nstorms) = 100. * real(islp)
             y_obs(nobs,nstorms) = lat_deg * deg2rad
          if ( GMT == 'GMT' ) then
!                                  Transfrom x from (-180 , 180) to (0, 360) then to radian
             if ( lon_deg < 0 ) then
                  x_obs(nobs,nstorms) = (360.+lon_deg) * deg2rad
             else
                  x_obs(nobs,nstorms) = (360.-lon_deg) * deg2rad
             endif
          elseif ( GMT == 'PAC' ) then   ! Pacific storms
             x_obs(nobs,nstorms) = lon_deg * deg2rad
          endif
     end select

  enddo

  endif

  close(unit)

  if(master) then
     write(*,*) 'TC vortex breeding: total storms=', nstorms
     if ( nstorms/=0 ) then
          do n=1,nstorms
             write(*, *) 'TC#',n, ' contains ',  nobs_tc(n),' observations'
          enddo
     endif
  endif

200  format(i3, 1x,f9.4, 1x, i2, 1x, i2, 1x, i2, 1x, a3, f5.1, 1x, f5.1, 1x, i3, 1x, i4, 1x, a19)

  end subroutine slp_obs_init

  !> The function 'calday' performs time interpolation:
  ! Julian day (0 to 365 for non-leap year)
  real function calday(year, month, day, hour, minute, sec)
! input:
    integer, intent(in):: year, month, day, hour
    integer, intent(in):: minute, sec
! Local:
      integer n, m, ds, nday
      real tsec
      integer days(12)
      data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

      ds = day - 1

      if( month /= 1 ) then
          do m=1, month-1
            if( m==2  .and. leap_year(year) ) then
                ds = ds + 29
            else
                ds = ds + days(m)
            endif
          enddo
      endif

      if ( leap_year(year_track_data) ) then
           nday = 366
      else
           nday = 365
      endif

      calday = real((year-year_track_data)*nday + ds)  + real(hour*3600 + minute*60 + sec)/86400.

  end function calday

!>@details The function 'leap_year' determines if year 'ny' is a leap year.
!>@author Shian-Jiann Lin
  logical function leap_year(ny)
  integer, intent(in):: ny
   integer ny00
!
! No leap years prior to 0000
!
      parameter ( ny00 = 0000 )   !< The threshold for starting leap-year

      if( ny >= ny00 ) then
         if( mod(ny,100) == 0. .and. mod(ny,400) == 0. ) then
             leap_year = .true.
         elseif( mod(ny,4) == 0. .and. mod(ny,100) /= 0.  ) then
             leap_year = .true.
         else
             leap_year = .false.
         endif
      else
          leap_year = .false.
      endif

  end function leap_year

 subroutine pmaxmin( qname, a, imax, jmax, fac )

      character(len=*)  qname
      integer imax, jmax
      integer i, j
      real a(imax,jmax)

      real qmin(jmax), qmax(jmax)
      real pmax, pmin
      real fac                     ! multiplication factor

      do j=1,jmax
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,imax
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
            pmax = qmax(1)
            pmin = qmin(1)
         do j=2,jmax
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

 end subroutine pmaxmin

!>@brief The subroutine 'del2_uv' filters the wind tendency.
 subroutine del2_uv(du, dv, cd, kmd, ntimes, bd, npx, npy, gridstruct, domain)
   integer, intent(in):: kmd
   integer, intent(in):: ntimes
   real,    intent(in):: cd            !< cd = K * da_min;   0 < K < 0.25
   type(fv_grid_bounds_type), intent(IN) :: bd
   real, intent(inout), dimension(is:ie,js:je,kmd):: du, dv
   integer, intent(IN) :: npx, npy
   type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: domain
! local:
  real(kind=R_GRID), pointer, dimension(:,:,:) :: vlon, vlat
   real, dimension(is:ie,js:je,kmd):: v1, v2, v3
   integer i,j,k

   vlon => gridstruct%vlon
   vlat => gridstruct%vlat

! transform to 3D Cartesian:
!$OMP parallel do default(none) shared(kmd,is,ie,js,je,v1,v2,v3,du,vlon,dv,vlat)
   do k=1,kmd
      do j=js,je
         do i=is,ie
            v1(i,j,k) = du(i,j,k)*vlon(i,j,1) + dv(i,j,k)*vlat(i,j,1)
            v2(i,j,k) = du(i,j,k)*vlon(i,j,2) + dv(i,j,k)*vlat(i,j,2)
            v3(i,j,k) = du(i,j,k)*vlon(i,j,3) + dv(i,j,k)*vlat(i,j,3)
         enddo
      enddo
   enddo

! Filter individual components as scalar:
   call del2_scalar( v1(is,js,1), cd, kmd, ntimes, bd, npx, npy, gridstruct, domain )
   call del2_scalar( v2(is,js,1), cd, kmd, ntimes, bd, npx, npy, gridstruct, domain )
   call del2_scalar( v3(is,js,1), cd, kmd, ntimes, bd, npx, npy, gridstruct, domain )

! Convert back to lat-lon components:
!$OMP parallel do default(none) shared(kmd,is,ie,js,je,du,dv,v1,v2,v3,vlon,vlat)
   do k=1,kmd
      do j=js,je
         do i=is,ie
            du(i,j,k) = v1(i,j,k)*vlon(i,j,1) + v2(i,j,k)*vlon(i,j,2) + v3(i,j,k)*vlon(i,j,3)
            dv(i,j,k) = v1(i,j,k)*vlat(i,j,1) + v2(i,j,k)*vlat(i,j,2) + v3(i,j,k)*vlat(i,j,3)
         enddo
      enddo
   enddo

 end subroutine del2_uv

!>@brief The subroutine 'del2_scalar' filters the physics tendency.
 subroutine del2_scalar(qdt, cd, kmd, nmax, bd,  npx, npy, gridstruct, domain)
   integer, intent(in):: kmd
   integer, intent(in):: nmax          !< must be no greater than 3
   real,    intent(in):: cd            !< cd = K * da_min;   0 < K < 0.25
   type(fv_grid_bounds_type), intent(IN) :: bd
   real, intent(inout):: qdt(is:ie,js:je,kmd)
   integer, intent(IN) :: npx, npy
   type(fv_grid_type), intent(IN), target :: gridstruct
  type(domain2d), intent(INOUT) :: domain
! local:
   real::  q(isd:ied,jsd:jed,kmd)
   real:: fx(isd:ied+1,jsd:jed), fy(isd:ied,jsd:jed+1)
   integer i,j,k, n, nt, ntimes
   real :: damp

  real, pointer, dimension(:,:) :: rarea, area
  real, pointer, dimension(:,:) :: sina_u, sina_v
  real, pointer, dimension(:,:,:) :: sin_sg

  real, pointer, dimension(:,:) :: dx, dy, rdxc, rdyc

  real(kind=R_GRID), pointer :: da_min

  logical, pointer :: bounded_domain, sw_corner, se_corner, nw_corner, ne_corner

   area => gridstruct%area
  rarea => gridstruct%rarea

  sina_u => gridstruct%sina_u
  sina_v => gridstruct%sina_v
  sin_sg => gridstruct%sin_sg

  dx     => gridstruct%dx
  dy     => gridstruct%dy
  rdxc   => gridstruct%rdxc
  rdyc   => gridstruct%rdyc

  da_min => gridstruct%da_min

  bounded_domain => gridstruct%bounded_domain
  sw_corner => gridstruct%sw_corner
  se_corner => gridstruct%se_corner
  nw_corner => gridstruct%nw_corner
  ne_corner => gridstruct%ne_corner

   ntimes = min(3, nmax)
   damp = cd * da_min
!$OMP parallel do default(none) shared(is,ie,js,je,kmd,q,qdt)
   do k=1,kmd
      ! mpp_update_domains fills neighbor halos, but physical outer halos of a
      ! bounded regional domain are untouched.  Initialize them before copying
      ! the compute region so extended filter stencils never read garbage.
      q(:,:,k) = 0.0
      do j=js,je
         do i=is,ie
            q(i,j,k) = qdt(i,j,k)
         enddo
      enddo
   enddo
                     call timing_on('COMM_TOTAL')
   call mpp_update_domains(q, domain, complete=.true.)
                     call timing_off('COMM_TOTAL')

   do n=1,ntimes

   nt = ntimes - n
!$OMP parallel do default(none) shared(is,ie,js,je,kmd,nt,dy,q,isd,jsd,npx,npy,bounded_domain,   &
!$OMP                                  bd,sw_corner,se_corner,nw_corner,ne_corner,       &
!$OMP                                  sina_u,rdxc,sin_sg,dx,rdyc,sina_v,qdt,damp,rarea, &
!$OMP                                  n,ntimes) &
!$OMP                          private(fx, fy)
   do k=1,kmd

      if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 1, bounded_domain, bd, &
           sw_corner, se_corner, nw_corner, ne_corner)
      do j=js-nt,je+nt
         do i=is-nt,ie+1+nt
            fx(i,j) = dy(i,j)*sina_u(i,j)*(q(i-1,j,k)-q(i,j,k))*rdxc(i,j)
         enddo
         if (is == 1) fx(is,j) = dy(is,j)*(q(is-1,j,k)-q(is,j,k))*rdxc(is,j)* &
            0.5*(sin_sg(1,j,1) + sin_sg(0,j,3))
         if (ie+1 == npx) fx(ie+1,j) = dy(ie+1,j)*(q(ie,j,k)-q(ie+1,j,k))*rdxc(ie+1,j)* &
            0.5*(sin_sg(npx,j,1) + sin_sg(npx-1,j,3))
      enddo

      if(nt>0) call copy_corners(q(isd,jsd,k), npx, npy, 2, bounded_domain, bd, &
           sw_corner, se_corner, nw_corner, ne_corner)

      do j=js-nt,je+1+nt
         if (j == 1 .OR. j == npy) then
            do i=is-nt,ie+nt
               fy(i,j) = dx(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j) &
                    *0.5*(sin_sg(i,j,2) + sin_sg(i,j-1,4) )
            enddo
         else
            do i=is-nt,ie+nt
               fy(i,j) = dx(i,j)*sina_v(i,j)*(q(i,j-1,k)-q(i,j,k))*rdyc(i,j)
            enddo
         end if
      enddo

      if ( nt==0 ) then
          do j=js,je
             do i=is,ie
                qdt(i,j,k) = q(i,j,k) + damp*rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
             enddo
          enddo
      else
          do j=js-nt,je+nt
             do i=is-nt,ie+nt
                q(i,j,k) = q(i,j,k) + damp*rarea(i,j)*(fx(i,j)-fx(i+1,j)+fy(i,j)-fy(i,j+1))
             enddo
          enddo
      endif
   enddo

   enddo
 end subroutine del2_scalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Adapted from GSI with simplification !!!!!!!!!!!!!!!!!!!!!
subroutine build_coarse_alpha_from_fine(alpha_x_fine, alpha_y_fine, &
                                        alpha_x_coarse, alpha_y_coarse, &
                                        isg, ieg, jsg, jeg, nxg_coarse, nyg_coarse, strideg)
  implicit none

  integer, intent(in) :: isg, ieg, jsg, jeg
  integer, intent(in) :: nxg_coarse, nyg_coarse, strideg

  real, intent(in)  :: alpha_x_fine(isg:ieg, jsg:jeg)
  real, intent(in)  :: alpha_y_fine(isg:ieg, jsg:jeg)

  real, intent(out) :: alpha_x_coarse(1:nxg_coarse, 1:nyg_coarse)
  real, intent(out) :: alpha_y_coarse(1:nxg_coarse, 1:nyg_coarse)

  integer :: ic, jc
  integer :: i_prev, i_curr, ii
  integer :: j_prev, j_curr, jj
  integer :: ifine, jfine
  real :: prod_alpha, aval

  alpha_x_coarse(:,:) = 0.0
  alpha_y_coarse(:,:) = 0.0

  ! ============================================================
  ! X direction
  !
  ! Convention:
  !   alpha_x_coarse(ic,jc) represents transition:
  !
  !       coarse ic-1  --->  coarse ic
  !
  !   Therefore alpha_x_coarse(2,jc) is the product of fine-grid
  !   alpha values from fine point corresponding to coarse 1
  !   to fine point corresponding to coarse 2.
  !
  !   alpha_x_coarse(1,jc) has no left neighbor. However, the
  !   current RF backward sweep uses alpha(1), so we set it equal
  !   to alpha(2) as a boundary convention.
  ! ============================================================

  do jc = 1, nyg_coarse

    jfine = jsg + (jc - 1) * strideg
    jfine = max(jsg, min(jfine, jeg))

    if (nxg_coarse == 1) then
      alpha_x_coarse(1,jc) = 0.0
    else

      do ic = 2, nxg_coarse

        i_prev = isg + (ic - 2) * strideg
        i_curr = isg + (ic - 1) * strideg

        i_prev = max(isg, min(i_prev, ieg))
        i_curr = max(isg, min(i_curr, ieg))

        prod_alpha = 1.0

        if (i_curr > i_prev) then
          do ii = i_prev + 1, i_curr
            aval = alpha_x_fine(ii, jfine)
            aval = min(0.999999, max(0.0, aval))
            prod_alpha = prod_alpha * aval
          enddo
        else
          aval = alpha_x_fine(i_curr, jfine)
          prod_alpha = min(0.999999, max(0.0, aval))
        endif

        alpha_x_coarse(ic,jc) = min(0.999999, max(0.0, prod_alpha))

      enddo

      ! Boundary convention for current backward sweep indexing.
      alpha_x_coarse(1,jc) = alpha_x_coarse(2,jc)

    endif

  enddo

  ! ============================================================
  ! Y direction
  !
  ! Convention:
  !   alpha_y_coarse(ic,jc) represents transition:
  !
  !       coarse jc-1  --->  coarse jc
  !
  !   alpha_y_coarse(ic,1) has no lower neighbor. Since the
  !   current RF backward sweep uses alpha(:,1), set it equal to
  !   alpha(:,2) as a boundary convention.
  ! ============================================================

  do ic = 1, nxg_coarse

    ifine = isg + (ic - 1) * strideg
    ifine = max(isg, min(ifine, ieg))

    if (nyg_coarse == 1) then
      alpha_y_coarse(ic,1) = 0.0
    else

      do jc = 2, nyg_coarse

        j_prev = jsg + (jc - 2) * strideg
        j_curr = jsg + (jc - 1) * strideg

        j_prev = max(jsg, min(j_prev, jeg))
        j_curr = max(jsg, min(j_curr, jeg))

        prod_alpha = 1.0

        if (j_curr > j_prev) then
          do jj = j_prev + 1, j_curr
            aval = alpha_y_fine(ifine, jj)
            aval = min(0.999999, max(0.0, aval))
            prod_alpha = prod_alpha * aval
          enddo
        else
          aval = alpha_y_fine(ifine, j_curr)
          prod_alpha = min(0.999999, max(0.0, aval))
        endif

        alpha_y_coarse(ic,jc) = min(0.999999, max(0.0, prod_alpha))

      enddo

      ! Boundary convention for current backward sweep indexing.
      alpha_y_coarse(ic,1) = alpha_y_coarse(ic,2)

    endif

  enddo

end subroutine build_coarse_alpha_from_fine

subroutine get_rf_domain_pelist(domain, pelist, npelist, root)
  use mpp_mod,         only: mpp_npes, mpp_pe, mpp_error, FATAL, &
                             mpp_send, mpp_recv, mpp_sync_self
  use mpp_domains_mod, only: domain2d, mpp_get_pelist
  implicit none

  type(domain2d), intent(in) :: domain
  integer, allocatable, intent(out) :: pelist(:)
  integer, intent(out) :: npelist, root

  integer :: npes_all
  integer :: n, m, pe, nvalid
  integer, allocatable :: pelist_tmp(:)
  integer, allocatable :: tile_of_pe(:)
  logical :: already_added
  logical, save :: printed_rf_pelist = .false.
  integer :: nshow, mype
  integer :: global_root

  mype = mpp_pe()
  npes_all = mpp_npes()
  global_root = 0

  !-----------------------------------------------------------------------
  ! For option 4, build the RF PE list from the actual tile number reported
  ! by each PE, rather than assuming contiguous PE blocks or relying on a
  ! domain2d pelist that may include all six cube faces.
  !
  !   - one tile  : all PEs with tile_num=1 are grouped together
  !   - six tiles : each tile gets the actual PEs whose nudge_tile_num matches
  !   - any N     : works as long as nudge_tile_num is assigned consistently
  !
  ! Option 1/3 is not affected because it does not enter RF.
  ! Option 2 falls through to the original/domain behavior below.
  !-----------------------------------------------------------------------
  if (nudge_option == 4) then

    allocate(tile_of_pe(0:npes_all-1))
    tile_of_pe(:) = -1

    if (mype == global_root) then
      tile_of_pe(mype) = nudge_tile_num

      do n = 0, npes_all-1
        if (n /= global_root) then
          call mpp_recv(tile_of_pe(n), glen=1, from_pe=n, block=.true.)
        endif
      enddo

      do n = 0, npes_all-1
        if (n /= global_root) then
          call mpp_send(tile_of_pe(0), plen=npes_all, to_pe=n)
        endif
      enddo

    else
      call mpp_send(nudge_tile_num, plen=1, to_pe=global_root)
      call mpp_recv(tile_of_pe(0), glen=npes_all, from_pe=global_root, block=.true.)
    endif

    call mpp_sync_self()

    if (nudge_tile_num <= 0) then
      call mpp_error(FATAL, 'get_rf_domain_pelist: option 4 requires valid nudge_tile_num on every PE')
    endif

    npelist = 0
    do n = 0, npes_all-1
      if (tile_of_pe(n) == nudge_tile_num) npelist = npelist + 1
    enddo

    if (npelist <= 0) then
      call mpp_error(FATAL, 'get_rf_domain_pelist: empty tile-based RF PE list')
    endif

    allocate(pelist(npelist))
    m = 0
    do n = 0, npes_all-1
      if (tile_of_pe(n) == nudge_tile_num) then
        m = m + 1
        pelist(m) = n
      endif
    enddo

    root = pelist(1)

    deallocate(tile_of_pe)

  else

    !--------------------------------------------------------------------
    ! Original/domain behavior for option 2 and other old RF paths.
    ! This intentionally avoids using the option-4 tile-number grouping.
    !--------------------------------------------------------------------
    allocate(pelist_tmp(npes_all))
    pelist_tmp(:) = -999999

    call mpp_get_pelist(domain, pelist_tmp)

    ! Count valid, unique PEs returned by the domain object.
    nvalid = 0
    do n = 1, npes_all
      pe = pelist_tmp(n)
      if (pe >= 0 .and. pe < npes_all) then
        already_added = .false.
        do m = 1, n - 1
          if (pelist_tmp(m) == pe) already_added = .true.
        enddo
        if (.not. already_added) nvalid = nvalid + 1
      endif
    enddo

    if (nvalid > 0) then
      npelist = nvalid
      allocate(pelist(npelist))

      nvalid = 0
      do n = 1, npes_all
        pe = pelist_tmp(n)
        if (pe >= 0 .and. pe < npes_all) then
          already_added = .false.
          do m = 1, nvalid
            if (pelist(m) == pe) already_added = .true.
          enddo
          if (.not. already_added) then
            nvalid = nvalid + 1
            pelist(nvalid) = pe
          endif
        endif
      enddo
    else
      ! Original fallback: all PEs participate in one RF group.
      npelist = npes_all
      allocate(pelist(npelist))
      do n = 1, npelist
        pelist(n) = n - 1
      enddo
    endif

    deallocate(pelist_tmp)

    if (npelist <= 0) then
      call mpp_error(FATAL, 'get_rf_domain_pelist: empty RF PE list')
    endif

    root = pelist(1)

  endif

end subroutine get_rf_domain_pelist

!=========================================================
! Helper: PE-list-local broadcast for RF work arrays.
! This avoids global collectives when different tiles use different roots.
!=========================================================
subroutine rf_pelist_broadcast_real(buf, nwords, root, pelist, npelist)
  use mpp_mod, only: mpp_pe, mpp_send, mpp_recv, mpp_sync_self
  implicit none

  integer, intent(in) :: nwords, root, npelist
  integer, intent(in) :: pelist(npelist)
  real,    intent(inout) :: buf(nwords)

  integer :: n, mype

  mype = mpp_pe()

  if (mype == root) then
    do n = 1, npelist
      if (pelist(n) /= root) then
        call mpp_send(buf(1), plen=nwords, to_pe=pelist(n))
      endif
    enddo
  else
    call mpp_recv(buf(1), glen=nwords, from_pe=root, block=.true.)
  endif

  call mpp_sync_self()

end subroutine rf_pelist_broadcast_real

!=========================================================
! Centralized cubed-sphere RF helpers for nudge_option=4.
!
! All MPI compute subdomains from all FV3 tiles are gathered to PE 0.
! PE 0 reconstructs the complete tile fields, uses the offline-verified
! fixed FV3 cube-edge connectivity, performs the RF on continuous periodic
! cube paths, and sends the filtered tile back to every PE.
!
! Edge numbering used below:
!   1 = WEST, 2 = EAST, 3 = SOUTH, 4 = NORTH
!=========================================================
subroutine rf_set_fixed_cube_topology(ntiles)
  implicit none
  integer, intent(in) :: ntiles

  rf_neighbor_tile(:,:) = 0
  rf_neighbor_edge(:,:) = 0
  rf_neighbor_reverse(:,:) = .false.

  ! This mapping is only defined for the standard six-tile global FV3 cube.
  ! Regional/nested configurations (including tile IDs 7/8) stay on the
  ! original option-1/2 paths and never use this topology.
  if (ntiles /= 6) return

  ! Verified offline from actual FV3 geolon/geolat geometry.
  ! Edge order: 1=WEST, 2=EAST, 3=SOUTH, 4=NORTH.
  !
  !              WEST        EAST        SOUTH       NORTH
  ! tile 1 :     5N(R)       2W          6N          3W(R)
  ! tile 2 :     1E          4S(R)       6E(R)       3S
  ! tile 3 :     1N(R)       4W          2N          5W(R)
  ! tile 4 :     3E          6S(R)       2E(R)       5S
  ! tile 5 :     3N(R)       6W          4N          1W(R)
  ! tile 6 :     5E          2S(R)       4E(R)       1S

  rf_neighbor_tile(1,:) = (/ 5, 2, 6, 3 /)
  rf_neighbor_tile(2,:) = (/ 1, 4, 6, 3 /)
  rf_neighbor_tile(3,:) = (/ 1, 4, 2, 5 /)
  rf_neighbor_tile(4,:) = (/ 3, 6, 2, 5 /)
  rf_neighbor_tile(5,:) = (/ 3, 6, 4, 1 /)
  rf_neighbor_tile(6,:) = (/ 5, 2, 4, 1 /)

  rf_neighbor_edge(1,:) = (/ 4, 1, 4, 1 /)
  rf_neighbor_edge(2,:) = (/ 2, 3, 2, 3 /)
  rf_neighbor_edge(3,:) = (/ 4, 1, 4, 1 /)
  rf_neighbor_edge(4,:) = (/ 2, 3, 2, 3 /)
  rf_neighbor_edge(5,:) = (/ 4, 1, 4, 1 /)
  rf_neighbor_edge(6,:) = (/ 2, 3, 2, 3 /)

  rf_neighbor_reverse(1,:) = (/ .true.,  .false., .false., .true.  /)
  rf_neighbor_reverse(2,:) = (/ .false., .true.,  .true.,  .false. /)
  rf_neighbor_reverse(3,:) = (/ .true.,  .false., .false., .true.  /)
  rf_neighbor_reverse(4,:) = (/ .false., .true.,  .true.,  .false. /)
  rf_neighbor_reverse(5,:) = (/ .true.,  .false., .false., .true.  /)
  rf_neighbor_reverse(6,:) = (/ .false., .true.,  .true.,  .false. /)
end subroutine rf_set_fixed_cube_topology

subroutine rf_build_symmetric_index_map(nfine, target_stride, ncoarse, cmap)
  use mpp_mod, only: mpp_error, FATAL
  implicit none
  integer, intent(in) :: nfine, target_stride
  integer, intent(out) :: ncoarse
  integer, allocatable, intent(out) :: cmap(:)

  integer :: ic, jc, stride_use, val
  real :: x

  if (nfine <= 1) then
    ncoarse = 1
    allocate(cmap(1))
    cmap(1) = 1
    return
  endif

  stride_use = max(1,target_stride)
  if (stride_use == 1) then
    ncoarse = nfine
    allocate(cmap(ncoarse))
    do ic=1,ncoarse
      cmap(ic)=ic
    enddo
    return
  endif

  ! Start from the requested spacing, then enforce the same parity as nfine.
  ! That allows exact reflection symmetry:
  !     cmap(ncoarse+1-ic) = nfine + 1 - cmap(ic)
  ! while preserving both endpoints.  For C384 and stride=4 this gives C98.
  ncoarse = (nfine-1 + stride_use-1)/stride_use + 1
  if (mod(ncoarse,2) /= mod(nfine,2)) ncoarse = min(nfine,ncoarse+1)
  ncoarse = max(2,min(nfine,ncoarse))

  allocate(cmap(ncoarse))
  cmap(:)=0

  do ic=1,ncoarse/2
    x = real(ic-1)*real(nfine-1)/real(ncoarse-1)
    val = 1 + nint(x)
    cmap(ic) = val
    jc = ncoarse + 1 - ic
    cmap(jc) = nfine + 1 - val
  enddo

  if (mod(ncoarse,2) == 1) cmap((ncoarse+1)/2)=(nfine+1)/2

  ! Defensive monotonicity check.
  do ic=2,ncoarse
    if (cmap(ic) <= cmap(ic-1)) then
      write(*,*) 'RF bad coarse index map:',nfine,target_stride,ncoarse,ic, &
                 cmap(ic-1),cmap(ic)
      call mpp_error(FATAL,'rf_build_symmetric_index_map: non-monotone map')
    endif
  enddo
end subroutine rf_build_symmetric_index_map

subroutine rf_build_interp_map(nfine, ncoarse, cmap, ilo, ihi, weight)
  implicit none
  integer, intent(in) :: nfine, ncoarse
  integer, intent(in) :: cmap(ncoarse)
  integer, intent(out) :: ilo(nfine), ihi(nfine)
  real, intent(out) :: weight(nfine)

  integer :: i, ic, den

  if (ncoarse <= 1) then
    ilo(:)=1
    ihi(:)=1
    weight(:)=0.0
    return
  endif

  ic=1
  do i=1,nfine
    do while (ic < ncoarse-1 .and. i > cmap(ic+1))
      ic=ic+1
    enddo

    if (i <= cmap(1)) then
      ilo(i)=1
      ihi(i)=1
      weight(i)=0.0
    else if (i >= cmap(ncoarse)) then
      ilo(i)=ncoarse
      ihi(i)=ncoarse
      weight(i)=0.0
    else
      ilo(i)=ic
      ihi(i)=ic+1
      den=cmap(ic+1)-cmap(ic)
      if (den > 0) then
        weight(i)=real(i-cmap(ic))/real(den)
      else
        weight(i)=0.0
      endif
      weight(i)=max(0.0,min(1.0,weight(i)))
    endif
  enddo
end subroutine rf_build_interp_map

subroutine build_coarse_alpha_from_map(alpha_x_fine, alpha_y_fine, &
                                       alpha_x_coarse, alpha_y_coarse, &
                                       nx, ny, nxc, nyc, imap, jmap)
  implicit none
  integer, intent(in) :: nx, ny, nxc, nyc
  integer, intent(in) :: imap(nxc), jmap(nyc)
  real, intent(in) :: alpha_x_fine(nx,ny), alpha_y_fine(nx,ny)
  real, intent(out) :: alpha_x_coarse(nxc,nyc), alpha_y_coarse(nxc,nyc)

  integer :: ic, jc, ii, jj, i_prev, i_curr, j_prev, j_curr
  integer :: ifine, jfine
  real :: prod_alpha, aval

  alpha_x_coarse(:,:)=0.0
  alpha_y_coarse(:,:)=0.0

  do jc=1,nyc
    jfine=jmap(jc)
    if (nxc == 1) then
      alpha_x_coarse(1,jc)=0.0
    else
      do ic=2,nxc
        i_prev=imap(ic-1)
        i_curr=imap(ic)
        prod_alpha=1.0
        do ii=i_prev+1,i_curr
          aval=min(0.999999,max(0.0,alpha_x_fine(ii,jfine)))
          prod_alpha=prod_alpha*aval
        enddo
        alpha_x_coarse(ic,jc)=min(0.999999,max(0.0,prod_alpha))
      enddo
      alpha_x_coarse(1,jc)=alpha_x_coarse(2,jc)
    endif
  enddo

  do ic=1,nxc
    ifine=imap(ic)
    if (nyc == 1) then
      alpha_y_coarse(ic,1)=0.0
    else
      do jc=2,nyc
        j_prev=jmap(jc-1)
        j_curr=jmap(jc)
        prod_alpha=1.0
        do jj=j_prev+1,j_curr
          aval=min(0.999999,max(0.0,alpha_y_fine(ifine,jj)))
          prod_alpha=prod_alpha*aval
        enddo
        alpha_y_coarse(ic,jc)=min(0.999999,max(0.0,prod_alpha))
      enddo
      alpha_y_coarse(ic,1)=alpha_y_coarse(ic,2)
    endif
  enddo
end subroutine build_coarse_alpha_from_map

subroutine rf_neighbor_ij(edge, depth, tang, nx, ny, ii, jj)
  implicit none
  integer, intent(in) :: edge, depth, tang, nx, ny
  integer, intent(out) :: ii, jj

  select case(edge)
  case(1)                    ! WEST: move inward with increasing i
    ii = depth;        jj = tang
  case(2)                    ! EAST: move inward with decreasing i
    ii = nx-depth+1;   jj = tang
  case(3)                    ! SOUTH: move inward with increasing j
    ii = tang;         jj = depth
  case(4)                    ! NORTH: move inward with decreasing j
    ii = tang;         jj = ny-depth+1
  case default
    ii = 1; jj = 1
  end select
end subroutine rf_neighbor_ij

subroutine rf_step_cube(nx, ny, ntiles, tile, i, j, dir, &
                        tile2, i2, j2, dir2)
  implicit none
  integer, intent(in) :: nx, ny, ntiles, tile, i, j, dir
  integer, intent(out) :: tile2, i2, j2, dir2

  integer :: nt, ne, tang, tang2, nedge

  tile2=tile
  i2=i
  j2=j
  dir2=dir

  select case(dir)
  case(1)                         ! WEST
    if (i > 1) then
      i2=i-1
      return
    endif
  case(2)                         ! EAST
    if (i < nx) then
      i2=i+1
      return
    endif
  case(3)                         ! SOUTH
    if (j > 1) then
      j2=j-1
      return
    endif
  case(4)                         ! NORTH
    if (j < ny) then
      j2=j+1
      return
    endif
  case default
    tile2=0
    return
  end select

  nt=rf_neighbor_tile(tile,dir)
  ne=rf_neighbor_edge(tile,dir)
  if (nt < 1 .or. nt > ntiles .or. ne < 1 .or. ne > 4) then
    tile2=0
    return
  endif

  if (dir <= 2) then
    tang=j
  else
    tang=i
  endif

  if (ne <= 2) then
    nedge=ny
  else
    nedge=nx
  endif

  tang2=tang
  if (rf_neighbor_reverse(tile,dir)) tang2=nedge+1-tang
  tang2=max(1,min(nedge,tang2))

  call rf_neighbor_ij(ne,1,tang2,nx,ny,i2,j2)
  tile2=nt

  select case(ne)
  case(1)                          ! enter WEST, continue EAST
    dir2=2
  case(2)                          ! enter EAST, continue WEST
    dir2=1
  case(3)                          ! enter SOUTH, continue NORTH
    dir2=4
  case(4)                          ! enter NORTH, continue SOUTH
    dir2=3
  end select
end subroutine rf_step_cube

subroutine rf_clear_cube_paths()
  implicit none

  if (allocated(rf_path_t)) deallocate(rf_path_t)
  if (allocated(rf_path_i)) deallocate(rf_path_i)
  if (allocated(rf_path_j)) deallocate(rf_path_j)
  if (allocated(rf_path_dir)) deallocate(rf_path_dir)
  if (allocated(rf_loop_start)) deallocate(rf_loop_start)
  if (allocated(rf_loop_len)) deallocate(rf_loop_len)
  if (allocated(rf_hit_count)) deallocate(rf_hit_count)

  rf_cube_paths_built=.false.
  rf_path_nx=0
  rf_path_ny=0
  rf_path_ntiles=0
  rf_path_nstate=0
  rf_path_nloop=0
  rf_path_maxloop=0
end subroutine rf_clear_cube_paths

subroutine rf_build_cube_paths(nx, ny, ntiles)
  use mpp_mod, only: mpp_error, FATAL
  implicit none
  integer, intent(in) :: nx, ny, ntiles

  integer :: maxstate, nstate, nloop, maxloop
  integer :: t, i, j, axis, dir, t2, i2, j2, dir2
  integer :: tnext, inext, jnext, axid, n
  integer :: start_t, start_i, start_j, start_axis
  logical, allocatable :: visited(:,:,:,:)

  if (rf_cube_paths_built) then
    if (rf_path_nx == nx .and. rf_path_ny == ny .and. &
        rf_path_ntiles == ntiles) return
  endif

  call rf_clear_cube_paths()

  if (ntiles /= 6) then
    call mpp_error(FATAL,'rf_build_cube_paths: continuous RF is only defined for six global FV3 tiles')
  endif
  if (nx /= ny) then
    call mpp_error(FATAL,'rf_build_cube_paths: six global FV3 faces must be square')
  endif

  call rf_set_fixed_cube_topology(ntiles)

  maxstate=2*nx*ny*ntiles
  allocate(visited(nx,ny,ntiles,2))
  allocate(rf_path_t(maxstate),rf_path_i(maxstate),rf_path_j(maxstate))
  allocate(rf_path_dir(maxstate))
  allocate(rf_loop_start(maxstate),rf_loop_len(maxstate))
  allocate(rf_hit_count(nx,ny,ntiles))

  visited=.false.
  rf_hit_count=0
  nstate=0
  nloop=0
  maxloop=0

  do t=1,ntiles
    do j=1,ny
      do i=1,nx
        do axis=1,2
          if (visited(i,j,t,axis)) cycle

          nloop=nloop+1
          rf_loop_start(nloop)=nstate+1
          start_t=t
          start_i=i
          start_j=j
          start_axis=axis
          if (axis == 1) then
            dir=2
          else
            dir=4
          endif
          t2=t
          i2=i
          j2=j
          n=0

          do
            if (dir <= 2) then
              axid=1
            else
              axid=2
            endif

            if (visited(i2,j2,t2,axid)) then
              if (t2 /= start_t .or. i2 /= start_i .or. j2 /= start_j .or. &
                  axid /= start_axis) then
                call mpp_error(FATAL,'rf_build_cube_paths: line merged into a previously visited line')
              endif
              exit
            endif

            visited(i2,j2,t2,axid)=.true.
            rf_hit_count(i2,j2,t2)=rf_hit_count(i2,j2,t2)+1

            nstate=nstate+1
            n=n+1
            if (nstate > maxstate) then
              call mpp_error(FATAL,'rf_build_cube_paths: state overflow')
            endif

            rf_path_t(nstate)=t2
            rf_path_i(nstate)=i2
            rf_path_j(nstate)=j2
            rf_path_dir(nstate)=dir

            call rf_step_cube(nx,ny,ntiles,t2,i2,j2,dir, &
                              tnext,inext,jnext,dir2)
            if (tnext < 1) then
              call mpp_error(FATAL,'rf_build_cube_paths: invalid cube topology step')
            endif
            t2=tnext
            i2=inext
            j2=jnext
            dir=dir2
          enddo

          rf_loop_len(nloop)=n
          maxloop=max(maxloop,n)
        enddo
      enddo
    enddo
  enddo

  if (nstate /= maxstate) then
    write(*,*) 'RF states expected/found=',maxstate,nstate
    call mpp_error(FATAL,'rf_build_cube_paths: not all coordinate states were traced')
  endif
  if (minval(rf_hit_count) /= 2 .or. maxval(rf_hit_count) /= 2) then
    write(*,*) 'RF hit-count min/max=',minval(rf_hit_count),maxval(rf_hit_count)
    call mpp_error(FATAL,'rf_build_cube_paths: every cell must belong to exactly two lines')
  endif

  rf_path_nx=nx
  rf_path_ny=ny
  rf_path_ntiles=ntiles
  rf_path_nstate=nstate
  rf_path_nloop=nloop
  rf_path_maxloop=maxloop
  rf_cube_paths_built=.true.

  deallocate(visited)

  if (nudge_debug) then
    write(*,'(a,3(i8,1x))') 'RF cached cube topology: loops states max_loop = ', &
         rf_path_nloop,rf_path_nstate,rf_path_maxloop
    flush(6)
  endif
end subroutine rf_build_cube_paths

subroutine rf_build_path_alpha(ax, ay, nx, ny, ntiles, path_alpha_out)
  implicit none
  integer, intent(in) :: nx, ny, ntiles
  real, intent(in) :: ax(nx,ny,ntiles), ay(nx,ny,ntiles)
  real, intent(out) :: path_alpha_out(rf_path_nstate)

  integer :: l, st, n, m, idx, idx2
  integer :: t, i, j, dir, t2, i2, j2, ne
  real :: a1, a2

  do l=1,rf_path_nloop
    st=rf_loop_start(l)
    n=rf_loop_len(l)
    do m=1,n
      idx=st+m-1
      if (m < n) then
        idx2=idx+1
      else
        idx2=st
      endif

      t=rf_path_t(idx)
      i=rf_path_i(idx)
      j=rf_path_j(idx)
      dir=rf_path_dir(idx)
      t2=rf_path_t(idx2)
      i2=rf_path_i(idx2)
      j2=rf_path_j(idx2)

      if (t2 == t) then
        select case(dir)
        case(1)
          a1=ax(i,j,t)
        case(2)
          a1=ax(i2,j2,t)
        case(3)
          a1=ay(i,j,t)
        case(4)
          a1=ay(i2,j2,t)
        end select
        path_alpha_out(idx)=min(0.999999,max(0.0,a1))
      else
        ne=rf_neighbor_edge(t,dir)
        if (dir <= 2) then
          a1=ax(i,j,t)
        else
          a1=ay(i,j,t)
        endif
        if (ne <= 2) then
          a2=ax(i2,j2,t2)
        else
          a2=ay(i2,j2,t2)
        endif
        a1=min(0.999999,max(0.0,a1))
        a2=min(0.999999,max(0.0,a2))
        path_alpha_out(idx)=sqrt(a1*a2)
      endif
    enddo
  enddo
end subroutine rf_build_path_alpha

subroutine rf1d_periodic_line(f, a, n, tmp, aa_line)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  integer, intent(in) :: n
  real, intent(inout) :: f(n)
  real, intent(in) :: a(n)
  real, intent(inout) :: tmp(:), aa_line(:)

  integer :: m
  real :: aa, bb, aprod, cterm, y0, z0

  if (n <= 1) return

  do m=1,n
    aa=a(m)
    if (ieee_is_nan(aa)) aa=0.0
    aa_line(m)=min(0.999999,max(0.0,aa))
  enddo

  ! Forward recursion on a closed line:
  !   y(m)=a(m)*y(m-1)+(1-a(m))*f(m), with y(0)=y(n).
  ! Solve the periodic closure analytically instead of choosing an endpoint.
  aprod=1.0
  cterm=0.0
  do m=1,n
    aa=aa_line(m)
    bb=1.0-aa
    cterm=aa*cterm+bb*f(m)
    aprod=aa*aprod
  enddo
  y0=cterm/max(1.0-aprod,1.0e-20)

  do m=1,n
    aa=aa_line(m)
    bb=1.0-aa
    tmp(m)=aa*y0+bb*f(m)
    y0=tmp(m)
  enddo

  ! Apply the exact adjoint of the variable-coefficient forward pass.
  ! If G is defined by
  !   y(m)=a(m)*y(m-1)+(1-a(m))*f(m),
  ! then G^T first solves
  !   z(m)=a(m+1)*z(m+1)+y(m)
  ! on the same periodic line, and only afterward multiplies by (1-a(m)).
  ! Folding (1-a(m+1)) into the reverse recurrence is correct only when
  ! alpha is constant and otherwise makes the filter non-symmetric.
  aprod=1.0
  cterm=0.0
  do m=n,1,-1
    if (m == n) then
      aa=aa_line(1)
    else
      aa=aa_line(m+1)
    endif
    cterm=aa*cterm+tmp(m)
    aprod=aa*aprod
  enddo
  z0=cterm/max(1.0-aprod,1.0e-20)

  do m=n,1,-1
    if (m == n) then
      aa=aa_line(1)
    else
      aa=aa_line(m+1)
    endif
    f(m)=aa*z0+tmp(m)
    z0=f(m)
  enddo
  do m=1,n
    f(m)=(1.0-aa_line(m))*f(m)
  enddo
end subroutine rf1d_periodic_line

subroutine rf_filter_cube_root(cube, ax, ay, nx, ny, nk, ntiles, nf)
  use mpp_mod, only: mpp_error, FATAL
  implicit none
  integer, intent(in) :: nx, ny, nk, ntiles, nf
  real, intent(inout) :: cube(nx,ny,nk,ntiles)
  real, intent(in) :: ax(nx,ny,ntiles), ay(nx,ny,ntiles)

  integer :: st, n, m, l, k, it, idx, t, i, j
  real, allocatable :: path_alpha_out(:)
  real, allocatable :: line(:), aline(:), tmp(:), aa_line(:)
  real, allocatable :: accum(:,:,:)

  if (ntiles /= 6) then
    call mpp_error(FATAL,'rf_filter_cube_root: continuous cube RF is only for six global FV3 tiles')
  endif
  if (nx /= ny) then
    call mpp_error(FATAL,'rf_filter_cube_root: six global FV3 faces must be square')
  endif

  ! Geometry is cached by active RF-grid dimensions.  rf_step_cube and
  ! rf_neighbor_ij are therefore used only when this coarse-grid topology
  ! is first built (or when nf_skip_grid changes the active dimensions).
  call rf_build_cube_paths(nx,ny,ntiles)

  allocate(path_alpha_out(rf_path_nstate))
  call rf_build_path_alpha(ax,ay,nx,ny,ntiles,path_alpha_out)

  allocate(line(rf_path_maxloop),aline(rf_path_maxloop))
  allocate(tmp(rf_path_maxloop),aa_line(rf_path_maxloop))
  allocate(accum(nx,ny,ntiles))

  do k=1,nk
    if (nudge_debug .and. (k == 1 .or. mod(k,10) == 0 .or. k == nk)) then
      write(*,'(a,i4,a,i4,a,2(es13.5,1x))') 'RF continuous cube level ',k,' / ',nk, &
           ' input min/max = ',minval(cube(:,:,k,:)),maxval(cube(:,:,k,:))
      flush(6)
    endif

    do it=1,nf
      accum(:,:,:)=0.0

      do l=1,rf_path_nloop
        st=rf_loop_start(l)
        n=rf_loop_len(l)

        do m=1,n
          idx=st+m-1
          line(m)=cube(rf_path_i(idx),rf_path_j(idx),k,rf_path_t(idx))
        enddo

        ! path_alpha_out(idx) is the transition idx -> idx+1.
        ! rf1d_periodic_line expects a(m) for m-1 -> m.
        aline(1)=path_alpha_out(st+n-1)
        do m=2,n
          aline(m)=path_alpha_out(st+m-2)
        enddo

        call rf1d_periodic_line(line(1:n),aline(1:n),n,tmp,aa_line)

        do m=1,n
          idx=st+m-1
          accum(rf_path_i(idx),rf_path_j(idx),rf_path_t(idx)) = &
               accum(rf_path_i(idx),rf_path_j(idx),rf_path_t(idx)) + line(m)
        enddo
      enddo

      do t=1,ntiles
        do j=1,ny
          do i=1,nx
            cube(i,j,k,t)=accum(i,j,t)/real(rf_hit_count(i,j,t))
          enddo
        enddo
      enddo
    enddo

    if (nudge_debug .and. (k == 1 .or. mod(k,10) == 0 .or. k == nk)) then
      write(*,'(a,i4,a,2(es13.5,1x))') 'RF continuous cube level ',k, &
           ' output min/max = ',minval(cube(:,:,k,:)),maxval(cube(:,:,k,:))
      flush(6)
    endif
  enddo

  deallocate(path_alpha_out,line,aline,tmp,aa_line,accum)
end subroutine rf_filter_cube_root

subroutine rf_debug_dump_slice(tag, field1, kmd, lon, lat, field2)
  use mpp_mod, only: mpp_pe, mpp_npes, mpp_send, mpp_recv, mpp_sync_self
  implicit none
  character(len=*), intent(in) :: tag
  integer, intent(in) :: kmd
  real, intent(in) :: field1(is:ie,js:je,1:kmd)
  real(kind=R_GRID), intent(in) :: lon(is:ie,js:je), lat(is:ie,js:je)
  real, intent(in), optional :: field2(is:ie,js:je,1:kmd)

  integer :: mype, npes, root, pe, p, i, j, idx, npt, nvar, nwords, pnwords, maxwords
  integer :: kk, iu
  integer :: local_meta(5)
  integer, allocatable :: meta(:,:)
  real, allocatable :: sbuf(:), rbuf(:)
  character(len=256) :: fname

  if (.not.rf_debug_dump) return

  kk = max(1,min(rf_debug_k,kmd))
  mype = mpp_pe(); npes = mpp_npes(); root = 0
  nvar = 3
  if (present(field2)) nvar = 4
  npt = (ie-is+1)*(je-js+1)
  nwords = nvar*npt
  allocate(sbuf(max(1,nwords)))

  idx=0
  do j=js,je
    do i=is,ie
      idx=idx+1; sbuf(idx)=real(lon(i,j))*180.0/pi
      idx=idx+1; sbuf(idx)=real(lat(i,j))*180.0/pi
      idx=idx+1; sbuf(idx)=field1(i,j,kk)
      if (present(field2)) then
        idx=idx+1; sbuf(idx)=field2(i,j,kk)
      endif
    enddo
  enddo

  local_meta=(/nudge_tile_num,is,ie,js,je/)
  if (mype == root) then
    allocate(meta(5,npes))
    meta(:,1)=local_meta
    do pe=1,npes-1
      call mpp_recv(meta(1,pe+1),glen=5,from_pe=pe,block=.true.)
    enddo
  else
    call mpp_send(local_meta(1),plen=5,to_pe=root)
  endif
  call mpp_sync_self()

  if (mype == root) then
    write(fname,'("sn_",a,"_k",i3.3,".csv")') trim(tag),kk
    open(newunit=iu,file=trim(fname),status='replace',action='write',form='formatted')
    if (present(field2)) then
      write(iu,'(a)') 'tile,i,j,lon_deg,lat_deg,value1,value2'
    else
      write(iu,'(a)') 'tile,i,j,lon_deg,lat_deg,value1'
    endif
    maxwords=1
    do p=1,npes
      pnwords=nvar*(meta(3,p)-meta(2,p)+1)*(meta(5,p)-meta(4,p)+1)
      maxwords=max(maxwords,pnwords)
    enddo
    allocate(rbuf(maxwords))

    do p=1,npes
      pnwords=nvar*(meta(3,p)-meta(2,p)+1)*(meta(5,p)-meta(4,p)+1)
      if (p == 1) then
        rbuf(1:pnwords)=sbuf(1:pnwords)
      else
        call mpp_recv(rbuf(1),glen=pnwords,from_pe=p-1,block=.true.)
      endif
      idx=0
      do j=meta(4,p),meta(5,p)
        do i=meta(2,p),meta(3,p)
          if (present(field2)) then
            write(iu,'(i0,a,i0,a,i0,a,4(es17.9,a))') meta(1,p),',',i,',',j,',', &
                 rbuf(idx+1),',',rbuf(idx+2),',',rbuf(idx+3),',',rbuf(idx+4),''
            idx=idx+4
          else
            write(iu,'(i0,a,i0,a,i0,a,3(es17.9,a))') meta(1,p),',',i,',',j,',', &
                 rbuf(idx+1),',',rbuf(idx+2),',',rbuf(idx+3),''
            idx=idx+3
          endif
        enddo
      enddo
    enddo
    close(iu)
    deallocate(rbuf,meta)
    write(*,'(a,a)') 'RF debug wrote ',trim(fname)
    flush(6)
  else
    call mpp_send(sbuf(1),plen=nwords,to_pe=root)
  endif
  call mpp_sync_self()
  deallocate(sbuf)
end subroutine rf_debug_dump_slice

subroutine build_rf_alpha_cube_from_L(gridstruct, Lx_m, Ly_m, bd, domain)
  use fv_arrays_mod,   only: fv_grid_type, fv_grid_bounds_type
  use mpp_mod,         only: mpp_pe, mpp_npes, mpp_send, mpp_recv, mpp_sync_self, mpp_error, FATAL
  use mpp_domains_mod, only: domain2d, mpp_get_global_domain, mpp_get_compute_domain
  implicit none

  type(domain2d),            intent(in) :: domain
  type(fv_grid_type),        intent(in), target :: gridstruct
  type(fv_grid_bounds_type), intent(in) :: bd
  real, intent(in) :: Lx_m, Ly_m

  integer :: mype, npes, root, isg, ieg, jsg, jeg, isc, iec, jsc, jec
  integer :: nx, ny, skipg, strideg, nxc, nyc
  integer :: p, pe, i, j, ii, jj, idx, nwords, maxwords
  integer :: n2, tile
  integer, allocatable :: meta(:,:), sizes(:), imap(:), jmap(:)
  integer :: local_meta(5)
  real, allocatable :: sendbuf(:), recvbuf(:)
  real, allocatable :: axfine(:,:,:), ayfine(:,:,:)
  real, pointer :: dx(:,:), dy(:,:)
  real :: axv, ayv
  real :: dx_min, dx_max, dy_min, dy_max
  real :: dx_sum, dy_sum, dx_sq_sum, dy_sq_sum
  integer :: cache_local_ok, cache_all_ok, peer_cache_ok

  mype=mpp_pe()
  npes=mpp_npes()
  root=0

  call mpp_get_global_domain(domain,isg,ieg,jsg,jeg)
  call mpp_get_compute_domain(domain,isc,iec,jsc,jec)
  dx=>gridstruct%dx
  dy=>gridstruct%dy

  skipg=max(0,nf_skip_grid)
  strideg=skipg+1
  nx=ieg-isg+1
  ny=jeg-jsg+1

  call rf_build_symmetric_index_map(nx,strideg,nxc,imap)
  call rf_build_symmetric_index_map(ny,strideg,nyc,jmap)

  dx_min=huge(1.0)
  dx_max=-huge(1.0)
  dy_min=huge(1.0)
  dy_max=-huge(1.0)
  dx_sum=0.0
  dy_sum=0.0
  dx_sq_sum=0.0
  dy_sq_sum=0.0
  do j=jsc,jec
    do i=isc,iec
      dx_min=min(dx_min,dx(i,j))
      dx_max=max(dx_max,dx(i,j))
      dy_min=min(dy_min,dy(i,j))
      dy_max=max(dy_max,dy(i,j))
      dx_sum=dx_sum+dx(i,j)
      dy_sum=dy_sum+dy(i,j)
      dx_sq_sum=dx_sq_sum+dx(i,j)*dx(i,j)
      dy_sq_sum=dy_sq_sum+dy(i,j)*dy(i,j)
    enddo
  enddo

  cache_local_ok=0
  if (rf_alpha_stored .and. rf_alpha_mode_stored == 2 .and. &
      rf_isg == isg .and. rf_ieg == ieg .and. rf_jsg == jsg .and. rf_jeg == jeg .and. &
      rf_alpha_stride == strideg .and. rf_alpha_tile_stored == nudge_tile_num .and. &
      rf_alpha_root_stored == root .and. &
      rf_alpha_isc_stored == isc .and. rf_alpha_iec_stored == iec .and. &
      rf_alpha_jsc_stored == jsc .and. rf_alpha_jec_stored == jec .and. &
      abs(rf_alpha_Lx_m_stored-Lx_m) <= 1.0e-6*max(1.0,abs(Lx_m)) .and. &
      abs(rf_alpha_Ly_m_stored-Ly_m) <= 1.0e-6*max(1.0,abs(Ly_m)) .and. &
      abs(rf_alpha_dx_min_stored-dx_min) <= 1.0e-6*max(1.0,abs(dx_min)) .and. &
      abs(rf_alpha_dx_max_stored-dx_max) <= 1.0e-6*max(1.0,abs(dx_max)) .and. &
      abs(rf_alpha_dy_min_stored-dy_min) <= 1.0e-6*max(1.0,abs(dy_min)) .and. &
      abs(rf_alpha_dy_max_stored-dy_max) <= 1.0e-6*max(1.0,abs(dy_max)) .and. &
      abs(rf_alpha_dx_sum_stored-dx_sum) <= 1.0e-6*max(1.0,abs(dx_sum)) .and. &
      abs(rf_alpha_dy_sum_stored-dy_sum) <= 1.0e-6*max(1.0,abs(dy_sum)) .and. &
      abs(rf_alpha_dx_sq_sum_stored-dx_sq_sum) <= 1.0e-6*max(1.0,abs(dx_sq_sum)) .and. &
      abs(rf_alpha_dy_sq_sum_stored-dy_sq_sum) <= 1.0e-6*max(1.0,abs(dy_sq_sum))) then
    cache_local_ok=1
  endif

  if (mype == root .and. cache_local_ok == 1) then
    if (.not.allocated(alpha_x_cube) .or. .not.allocated(alpha_y_cube)) then
      cache_local_ok=0
    else if (size(alpha_x_cube,1) /= nxc .or. size(alpha_x_cube,2) /= nyc .or. &
             size(alpha_x_cube,3) /= 6 .or. size(alpha_y_cube,1) /= nxc .or. &
             size(alpha_y_cube,2) /= nyc .or. size(alpha_y_cube,3) /= 6 .or. &
             rf_cube_ntiles /= 6) then
      cache_local_ok=0
    endif
  endif

  ! Alpha lives on PE 0; require all global PEs to agree before reusing it.
  if (mype == root) then
    cache_all_ok=cache_local_ok
    do pe=1,npes-1
      call mpp_recv(peer_cache_ok,glen=1,from_pe=pe,block=.true.)
      cache_all_ok=min(cache_all_ok,peer_cache_ok)
    enddo
    do pe=1,npes-1
      call mpp_send(cache_all_ok,plen=1,to_pe=pe)
    enddo
  else
    call mpp_send(cache_local_ok,plen=1,to_pe=root)
    call mpp_recv(cache_all_ok,glen=1,from_pe=root,block=.true.)
  endif
  call mpp_sync_self()

  if (cache_all_ok == 1) then
    deallocate(imap,jmap)
    return
  endif

  ! Alpha is rebuilt infrequently (only when L/grid/stride changes), so keep
  ! the robust fine-grid gather here.  The large 3-D field gather below is
  ! sampled BEFORE MPI communication, which is where skip-grid saves the
  ! dominant memory and communication cost.
  n2=(iec-isc+1)*(jec-jsc+1)
  nwords=2*n2
  allocate(sendbuf(max(1,nwords)))

  idx=0
  do j=jsc,jec
    do i=isc,iec
      idx=idx+1
      axv=exp(-dx(i,j)/max(Lx_m,1.0))
      sendbuf(idx)=min(0.999999,max(0.0,axv))
    enddo
  enddo
  do j=jsc,jec
    do i=isc,iec
      idx=idx+1
      ayv=exp(-dy(i,j)/max(Ly_m,1.0))
      sendbuf(idx)=min(0.999999,max(0.0,ayv))
    enddo
  enddo

  local_meta=(/nudge_tile_num,isc,iec,jsc,jec/)
  if (mype == root) then
    allocate(meta(5,npes),sizes(npes))
    meta(:,1)=local_meta
    sizes(1)=nwords
    do pe=1,npes-1
      call mpp_recv(meta(1,pe+1),glen=5,from_pe=pe,block=.true.)
      call mpp_recv(sizes(pe+1),glen=1,from_pe=pe,block=.true.)
    enddo
  else
    call mpp_send(local_meta(1),plen=5,to_pe=root)
    call mpp_send(nwords,plen=1,to_pe=root)
  endif
  call mpp_sync_self()

  if (mype == root) then
    ! This routine is entered only for nudge_option=4 with the continuous
    ! global RF enabled.  Require the standard global tiles 1..6 here.
    do p=1,npes
      tile=meta(1,p)
      if (tile < 1 .or. tile > 6) then
        write(*,*) 'RF unexpected global tile id=',tile,' pe=',p-1
        call mpp_error(FATAL,'build_rf_alpha_cube_from_L: continuous RF requires tiles 1..6')
      endif
    enddo
    do tile=1,6
      if (count(meta(1,:) == tile) == 0) then
        write(*,*) 'RF missing global tile=',tile
        call mpp_error(FATAL,'build_rf_alpha_cube_from_L: missing one of six global tiles')
      endif
    enddo
    rf_cube_ntiles=6

    allocate(axfine(nx,ny,6),ayfine(nx,ny,6))
    axfine=0.0
    ayfine=0.0
    maxwords=max(1,maxval(sizes))
    allocate(recvbuf(maxwords))

    do p=1,npes
      if (p == 1) then
        recvbuf(1:sizes(p))=sendbuf(1:sizes(p))
      else
        call mpp_recv(recvbuf(1),glen=sizes(p),from_pe=p-1,block=.true.)
      endif

      n2=(meta(3,p)-meta(2,p)+1)*(meta(5,p)-meta(4,p)+1)
      idx=0
      do j=meta(4,p),meta(5,p)
        do i=meta(2,p),meta(3,p)
          idx=idx+1
          ii=i-isg+1
          jj=j-jsg+1
          axfine(ii,jj,meta(1,p))=recvbuf(idx)
        enddo
      enddo
      do j=meta(4,p),meta(5,p)
        do i=meta(2,p),meta(3,p)
          idx=idx+1
          ii=i-isg+1
          jj=j-jsg+1
          ayfine(ii,jj,meta(1,p))=recvbuf(idx)
        enddo
      enddo
    enddo
    deallocate(recvbuf)

  else
    call mpp_send(sendbuf(1),plen=nwords,to_pe=root)
  endif
  call mpp_sync_self()
  deallocate(sendbuf)

  if (mype == root) then
    call rf_set_fixed_cube_topology(6)

    if (allocated(alpha_x_cube)) deallocate(alpha_x_cube)
    if (allocated(alpha_y_cube)) deallocate(alpha_y_cube)
    allocate(alpha_x_cube(nxc,nyc,6),alpha_y_cube(nxc,nyc,6))

    if (strideg == 1) then
      alpha_x_cube(:,:,:)=axfine(:,:,:)
      alpha_y_cube(:,:,:)=ayfine(:,:,:)
    else
      do tile=1,6
        call build_coarse_alpha_from_map(axfine(:,:,tile),ayfine(:,:,tile), &
             alpha_x_cube(:,:,tile),alpha_y_cube(:,:,tile), &
             nx,ny,nxc,nyc,imap,jmap)
      enddo
    endif

    if (nudge_debug) write(*,'(a,4(i6,1x))') &
         'RF centralized alpha nx ny target_stride tiles = ',nxc,nyc,strideg,6
    if (nudge_debug .and. strideg > 1) then
      write(*,'(a,4(i6,1x))') 'RF fine/coarse X endpoints = ', &
           imap(1),imap(nxc),nxc,nx
      write(*,'(a,4(i6,1x))') 'RF fine/coarse Y endpoints = ', &
           jmap(1),jmap(nyc),nyc,ny
    endif

    deallocate(axfine,ayfine,meta,sizes)
  endif

  deallocate(imap,jmap)

  rf_isg=isg
  rf_ieg=ieg
  rf_jsg=jsg
  rf_jeg=jeg
  rf_alpha_stride=strideg
  rf_alpha_Lx_m_stored=Lx_m
  rf_alpha_Ly_m_stored=Ly_m
  rf_alpha_mode_stored=2
  rf_alpha_tile_stored=nudge_tile_num
  rf_alpha_root_stored=root
  rf_alpha_isc_stored=isc
  rf_alpha_iec_stored=iec
  rf_alpha_jsc_stored=jsc
  rf_alpha_jec_stored=jec
  rf_alpha_dx_min_stored=dx_min
  rf_alpha_dx_max_stored=dx_max
  rf_alpha_dy_min_stored=dy_min
  rf_alpha_dy_max_stored=dy_max
  rf_alpha_dx_sum_stored=dx_sum
  rf_alpha_dy_sum_stored=dy_sum
  rf_alpha_dx_sq_sum_stored=dx_sq_sum
  rf_alpha_dy_sq_sum_stored=dy_sq_sum
  rf_alpha_stored=.true.
end subroutine build_rf_alpha_cube_from_L

subroutine rf2d_apply_cube_global(field_inout, bd, nf, isg, ieg, jsg, jeg, kmd, domain)
  use mpp_mod, only: mpp_pe, mpp_npes, mpp_send, mpp_recv, mpp_sync_self, mpp_error, FATAL
  use mpp_domains_mod, only: domain2d
  use fv_arrays_mod, only: fv_grid_bounds_type
  implicit none

  type(domain2d), intent(inout) :: domain
  type(fv_grid_bounds_type), intent(in) :: bd
  integer, intent(in) :: nf,isg,ieg,jsg,jeg,kmd
  real, intent(inout) :: field_inout(is:ie,js:je,1:kmd)

  integer :: mype,npes,root,skipg,strideg,nx,ny,nxc,nyc,nact_loc
  integer :: i,j,k,p,pe,idx,nwords,maxwords,tile,ic,jc
  integer :: cis,cie,cjs,cje,ni_s,nj_s,ig,jg
  integer :: fis,fie,fjs,fje,ric1,ric2,rjc1,rjc2
  integer :: ni_r,nj_r,nrecv,maxscatter,sp,pp
  integer :: pc1,pc2,pj1,pj2,pni,pnj,pwords
  integer, allocatable :: meta(:,:),sizes(:),imap(:),jmap(:)
  integer, allocatable :: coarse_hit(:,:,:)
  integer, allocatable :: ilo_i(:),ihi_i(:),ilo_j(:),ihi_j(:)
  integer :: local_meta(9)
  real, allocatable :: wi_i(:),wj_j(:)
  real, allocatable :: sendbuf(:),recvbuf(:),cube(:,:,:,:)
  real, allocatable :: face(:,:,:),scatterbuf(:)
  real :: wi,wj

  mype=mpp_pe()
  npes=mpp_npes()
  root=0
  skipg=max(0,nf_skip_grid)
  strideg=skipg+1
  nx=ieg-isg+1
  ny=jeg-jsg+1
  nact_loc=min(nact,kmd)

  call rf_build_symmetric_index_map(nx,strideg,nxc,imap)
  call rf_build_symmetric_index_map(ny,strideg,nyc,jmap)

  allocate(ilo_i(nx),ihi_i(nx),wi_i(nx))
  allocate(ilo_j(ny),ihi_j(ny),wj_j(ny))
  call rf_build_interp_map(nx,nxc,imap,ilo_i,ihi_i,wi_i)
  call rf_build_interp_map(ny,nyc,jmap,ilo_j,ihi_j,wj_j)

  if (mype == root) then
    if (.not.allocated(alpha_x_cube) .or. .not.allocated(alpha_y_cube)) then
      call mpp_error(FATAL,'rf2d_apply_cube_global: centralized alpha cube is not built')
    endif
    if (size(alpha_x_cube,1) /= nxc .or. size(alpha_x_cube,2) /= nyc) then
      call mpp_error(FATAL,'rf2d_apply_cube_global: alpha/field coarse-grid mismatch')
    endif
    if (rf_cube_ntiles /= 6) then
      call mpp_error(FATAL,'rf2d_apply_cube_global: continuous RF requires six global tiles')
    endif
  endif

  ! Fine-grid local compute bounds expressed in 1..nx / 1..ny face indices.
  fis=is-isg+1
  fie=ie-isg+1
  fjs=js-jsg+1
  fje=je-jsg+1

  ! Identify coarse map points owned by this PE.  The 3-D field is sampled
  ! here, BEFORE MPI communication to PE 0.
  cis=1
  cie=0
  do ic=1,nxc
    ig=isg+imap(ic)-1
    if (ig >= is .and. ig <= ie) then
      if (cie < cis) cis=ic
      cie=ic
    endif
  enddo

  cjs=1
  cje=0
  do jc=1,nyc
    jg=jsg+jmap(jc)-1
    if (jg >= js .and. jg <= je) then
      if (cje < cjs) cjs=jc
      cje=jc
    endif
  enddo

  if (cie >= cis) then
    ni_s=cie-cis+1
  else
    ni_s=0
  endif
  if (cje >= cjs) then
    nj_s=cje-cjs+1
  else
    nj_s=0
  endif

  if (ni_s > 0 .and. nj_s > 0) then
    nwords=ni_s*nj_s*nact_loc
  else
    nwords=0
  endif

  allocate(sendbuf(max(1,nwords)))
  idx=0
  if (nwords > 0) then
    do k=1,nact_loc
      do jc=cjs,cje
        jg=jsg+jmap(jc)-1
        do ic=cis,cie
          ig=isg+imap(ic)-1
          idx=idx+1
          sendbuf(idx)=field_inout(ig,jg,klist(k))
        enddo
      enddo
    enddo
  endif

  ! Metadata:
  ! 1 tile
  ! 2:5 local fine bounds on the complete face
  ! 6:9 locally owned coarse-map rectangle
  local_meta=(/nudge_tile_num,fis,fie,fjs,fje,cis,cie,cjs,cje/)

  if (mype == root) then
    allocate(meta(9,npes),sizes(npes))
    meta(:,1)=local_meta
    sizes(1)=nwords
    do pe=1,npes-1
      call mpp_recv(meta(1,pe+1),glen=9,from_pe=pe,block=.true.)
      call mpp_recv(sizes(pe+1),glen=1,from_pe=pe,block=.true.)
    enddo
  else
    call mpp_send(local_meta(1),plen=9,to_pe=root)
    call mpp_send(nwords,plen=1,to_pe=root)
  endif
  call mpp_sync_self()

  if (mype == root) then
    do p=1,npes
      if (meta(1,p) < 1 .or. meta(1,p) > 6) then
        call mpp_error(FATAL,'rf2d_apply_cube_global: unexpected tile id in six-tile global gather')
      endif
    enddo

    allocate(cube(nxc,nyc,nact_loc,6))
    allocate(coarse_hit(nxc,nyc,6))
    cube=0.0
    coarse_hit=0
    maxwords=max(1,maxval(sizes))
    allocate(recvbuf(maxwords))

    do p=1,npes
      if (sizes(p) > 0) then
        if (p == 1) then
          recvbuf(1:sizes(p))=sendbuf(1:sizes(p))
        else
          call mpp_recv(recvbuf(1),glen=sizes(p),from_pe=p-1,block=.true.)
        endif

        idx=0
        do k=1,nact_loc
          do jc=meta(8,p),meta(9,p)
            do ic=meta(6,p),meta(7,p)
              idx=idx+1
              cube(ic,jc,k,meta(1,p))=recvbuf(idx)
              if (k == 1) coarse_hit(ic,jc,meta(1,p)) = coarse_hit(ic,jc,meta(1,p)) + 1
            enddo
          enddo
        enddo
      endif
    enddo
    deallocate(recvbuf)

    ! Every retained coarse-grid point on every tile must have been supplied
    ! by exactly one PE.  This validates the global symmetric skip-grid gather
    ! independently of the MPI decomposition.
    if (minval(coarse_hit) /= 1 .or. maxval(coarse_hit) /= 1) then
      write(*,'(a,2(i8,1x))') 'RF coarse gather hit min/max = ', &
           minval(coarse_hit),maxval(coarse_hit)
      call mpp_error(FATAL,'rf2d_apply_cube_global: coarse gather has missing or duplicated points')
    endif
    if (nudge_debug) then
      write(*,'(a,2(i8,1x))') 'RF coarse gather hit min/max = ', &
           minval(coarse_hit),maxval(coarse_hit)
      flush(6)
    endif
  else
    if (nwords > 0) call mpp_send(sendbuf(1),plen=nwords,to_pe=root)
  endif
  call mpp_sync_self()
  deallocate(sendbuf)

  if (mype == root) then
    if (nudge_debug) write(*,'(a,4(i6,1x))') &
         'RF centralized cube nx ny levels target_stride = ',nxc,nyc,nact_loc,strideg
    call rf_filter_cube_root(cube,alpha_x_cube,alpha_y_cube,nxc,nyc,nact_loc,6,nf)
    if (nudge_debug) then
      write(*,'(a)') 'RF centralized cube filter COMPLETE'
      flush(6)
    endif
  endif

  ! Each PE receives only the coarse rectangle needed to interpolate its
  ! own fine-grid compute domain.  This replaces the previous whole-face
  ! return to every PE.
  ric1=ilo_i(fis)
  ric2=ihi_i(fie)
  rjc1=ilo_j(fjs)
  rjc2=ihi_j(fje)
  ni_r=ric2-ric1+1
  nj_r=rjc2-rjc1+1
  nrecv=ni_r*nj_r*nact_loc
  allocate(face(ric1:ric2,rjc1:rjc2,nact_loc))

  if (mype == root) then
    maxscatter=max(1,nrecv)
    do p=1,npes
      pc1=ilo_i(meta(2,p))
      pc2=ihi_i(meta(3,p))
      pj1=ilo_j(meta(4,p))
      pj2=ihi_j(meta(5,p))
      pni=pc2-pc1+1
      pnj=pj2-pj1+1
      maxscatter=max(maxscatter,pni*pnj*nact_loc)
    enddo
    allocate(scatterbuf(maxscatter))

    do p=1,npes
      tile=meta(1,p)
      pc1=ilo_i(meta(2,p))
      pc2=ihi_i(meta(3,p))
      pj1=ilo_j(meta(4,p))
      pj2=ihi_j(meta(5,p))
      pni=pc2-pc1+1
      pnj=pj2-pj1+1
      pwords=pni*pnj*nact_loc

      idx=0
      do k=1,nact_loc
        do jc=pj1,pj2
          do ic=pc1,pc2
            idx=idx+1
            scatterbuf(idx)=cube(ic,jc,k,tile)
          enddo
        enddo
      enddo

      if (p == 1) then
        idx=0
        do k=1,nact_loc
          do jc=rjc1,rjc2
            do ic=ric1,ric2
              idx=idx+1
              face(ic,jc,k)=scatterbuf(idx)
            enddo
          enddo
        enddo
      else
        pe=p-1
        if (nudge_debug .and. (pe == 1 .or. mod(pe,50) == 0 .or. pe == npes-1)) then
          write(*,'(a,i6,a,i6,a,i2,a,i8)') 'RF scatter subdomain to PE ',pe, &
               ' / ',npes-1,' tile=',tile,' words=',pwords
          flush(6)
        endif
        call mpp_send(scatterbuf(1),plen=pwords,to_pe=pe)
        ! mpp_send may not have consumed scatterbuf on return.  Because the
        ! next PE gets a DIFFERENT packed rectangle, complete this send before
        ! scatterbuf is overwritten for the next PE.
        call mpp_sync_self()
      endif
    enddo

    deallocate(scatterbuf)
    if (nudge_debug) then
      write(*,'(a)') 'RF centralized subdomain scatter COMPLETE'
      flush(6)
    endif
  else
    call mpp_recv(face(ric1,rjc1,1),glen=nrecv,from_pe=root,block=.true.)
  endif
  call mpp_sync_self()

  ! Interpolate the endpoint-preserving coarse map back to the local fine
  ! compute domain.  With stride=1 this reduces exactly to direct copy.
  do k=1,nact_loc
    do j=js,je
      jc=ilo_j(j-jsg+1)
      pnj=ihi_j(j-jsg+1)
      wj=wj_j(j-jsg+1)
      do i=is,ie
        ic=ilo_i(i-isg+1)
        pni=ihi_i(i-isg+1)
        wi=wi_i(i-isg+1)

        field_inout(i,j,klist(k)) = &
          (1.0-wi)*(1.0-wj)*face(ic ,jc ,k) + &
               wi *(1.0-wj)*face(pni,jc ,k) + &
          (1.0-wi)*     wj *face(ic ,pnj,k) + &
               wi *     wj *face(pni,pnj,k)
      enddo
    enddo
  enddo

  deallocate(face,imap,jmap,ilo_i,ihi_i,ilo_j,ihi_j,wi_i,wj_j)
  if (mype == root) deallocate(cube,meta,sizes,coarse_hit)
end subroutine rf2d_apply_cube_global

subroutine build_rf_alpha_from_L(gridstruct, Lx_m, Ly_m, bd, domain)
  use fv_arrays_mod,   only: fv_grid_type, fv_grid_bounds_type
  use mpp_mod,         only: mpp_pe, mpp_send, mpp_recv, mpp_sync_self, mpp_npes
  use mpp_domains_mod, only: domain2d, mpp_get_global_domain, mpp_get_compute_domain, mpp_get_pelist

  implicit none

  type(domain2d),            intent(in) :: domain
  type(fv_grid_type),        intent(in), target :: gridstruct
  type(fv_grid_bounds_type), intent(in) :: bd
  real,                      intent(in) :: Lx_m, Ly_m

  ! local alpha fields (only compute part will be used)
  real :: alpha_x(is:ie, js:je)
  real :: alpha_y(is:ie, js:je)

  real, pointer :: dx(:,:), dy(:,:)
  integer :: i, j

  integer :: isg, ieg, jsg, jeg
  integer :: isc, iec, jsc, jec
  integer :: ni, nj, n2, send_count
  integer :: mype, root, pe_root

  integer :: npelist, n
  integer, allocatable :: pelist(:)
  integer :: Ldims(4)
  integer, allocatable :: Gdims(:)     ! 4*npelist
  integer, allocatable :: LsizeS(:)    ! per-PE packed size
  integer, allocatable :: displ(:)     ! displacement into recvbuf

  real, allocatable :: sendbuf(:)
  real, allocatable :: recvbuf(:)      ! only meaningful on root
  integer :: tot

  ! Active RF alpha grid.  If nf_skip_grid > 0, alpha_x_glob/alpha_y_glob
  !     are stored on the coarse grid directly, not as fine-grid arrays.
  integer :: skipg, strideg
  integer :: nxg_alpha, nyg_alpha
  real, allocatable :: alpha_x_fine_root(:,:), alpha_y_fine_root(:,:)

  integer :: cache_local_ok, cache_all_ok, peer_cache_ok
  real :: dx_min, dx_max, dy_min, dy_max
  real :: dx_sum, dy_sum, dx_sq_sum, dy_sq_sum

  real :: ax, ay
  integer :: p, i1, i2, j1, j2, nij, off, kk
  integer :: base


  ! -----------------------------
  ! Domain / PE info
  ! -----------------------------
  ! Comparison switch for option 4:
  !   rf_cube_continuous=.true.  -> centralized six-tile continuous RF.
  !   rf_cube_continuous=.false. -> ORIGINAL independent-tile RF below;
  !                                 each tile is filtered independently.
  ! The CSV dumps are outside this branch, so before/after diagnostics are
  ! produced identically for either mode.
  if (nudge_option == 4 .and. rf_cube_continuous) then
    call build_rf_alpha_cube_from_L(gridstruct, Lx_m, Ly_m, bd, domain)
    return
  endif

  mype = mpp_pe()

  call get_rf_domain_pelist(domain, pelist, npelist, root)

  dx => gridstruct%dx
  dy => gridstruct%dy

  call mpp_get_global_domain (domain, isg, ieg, jsg, jeg)
  call mpp_get_compute_domain(domain, isc, iec, jsc, jec)

  skipg = max(0, nf_skip_grid)
  strideg = skipg + 1

  if (skipg == 0) then
    nxg_alpha = ieg - isg + 1
    nyg_alpha = jeg - jsg + 1
  else
    nxg_alpha = (ieg - isg) / strideg + 1
    nyg_alpha = (jeg - jsg) / strideg + 1
  endif

  dx_min = huge(1.0)
  dx_max = -huge(1.0)
  dy_min = huge(1.0)
  dy_max = -huge(1.0)
  dx_sum = 0.0
  dy_sum = 0.0
  dx_sq_sum = 0.0
  dy_sq_sum = 0.0
  do j = jsc, jec
    do i = isc, iec
      dx_min = min(dx_min, dx(i,j))
      dx_max = max(dx_max, dx(i,j))
      dy_min = min(dy_min, dy(i,j))
      dy_max = max(dy_max, dy(i,j))
      dx_sum = dx_sum + dx(i,j)
      dy_sum = dy_sum + dy(i,j)
      dx_sq_sum = dx_sq_sum + dx(i,j) * dx(i,j)
      dy_sq_sum = dy_sq_sum + dy(i,j) * dy(i,j)
    enddo
  enddo

  cache_local_ok = 0
  if (rf_alpha_stored .and. rf_alpha_mode_stored == 1 .and. &
      rf_isg == isg .and. rf_ieg == ieg .and. &
      rf_jsg == jsg .and. rf_jeg == jeg .and. &
      rf_alpha_stride == strideg .and. &
      rf_alpha_tile_stored == nudge_tile_num .and. &
      rf_alpha_root_stored == root .and. &
      rf_alpha_isc_stored == isc .and. rf_alpha_iec_stored == iec .and. &
      rf_alpha_jsc_stored == jsc .and. rf_alpha_jec_stored == jec .and. &
      abs(rf_alpha_Lx_m_stored - Lx_m) <= 1.0e-6 * max(1.0,abs(Lx_m)) .and. &
      abs(rf_alpha_Ly_m_stored - Ly_m) <= 1.0e-6 * max(1.0,abs(Ly_m)) .and. &
      abs(rf_alpha_dx_min_stored - dx_min) <= 1.0e-6 * max(1.0,abs(dx_min)) .and. &
      abs(rf_alpha_dx_max_stored - dx_max) <= 1.0e-6 * max(1.0,abs(dx_max)) .and. &
      abs(rf_alpha_dy_min_stored - dy_min) <= 1.0e-6 * max(1.0,abs(dy_min)) .and. &
      abs(rf_alpha_dy_max_stored - dy_max) <= 1.0e-6 * max(1.0,abs(dy_max)) .and. &
      abs(rf_alpha_dx_sum_stored - dx_sum) <= 1.0e-6 * max(1.0,abs(dx_sum)) .and. &
      abs(rf_alpha_dy_sum_stored - dy_sum) <= 1.0e-6 * max(1.0,abs(dy_sum)) .and. &
      abs(rf_alpha_dx_sq_sum_stored - dx_sq_sum) <= 1.0e-6 * max(1.0,abs(dx_sq_sum)) .and. &
      abs(rf_alpha_dy_sq_sum_stored - dy_sq_sum) <= 1.0e-6 * max(1.0,abs(dy_sq_sum))) then
    cache_local_ok = 1
  endif

  if (mype == root .and. cache_local_ok == 1) then
    if (.not.allocated(alpha_x_glob) .or. .not.allocated(alpha_y_glob)) then
      cache_local_ok = 0
    else if (skipg == 0) then
      if (lbound(alpha_x_glob,1) /= isg .or. ubound(alpha_x_glob,1) /= ieg .or. &
          lbound(alpha_x_glob,2) /= jsg .or. ubound(alpha_x_glob,2) /= jeg .or. &
          lbound(alpha_y_glob,1) /= isg .or. ubound(alpha_y_glob,1) /= ieg .or. &
          lbound(alpha_y_glob,2) /= jsg .or. ubound(alpha_y_glob,2) /= jeg) then
        cache_local_ok = 0
      endif
    else
      if (size(alpha_x_glob,1) /= nxg_alpha .or. size(alpha_x_glob,2) /= nyg_alpha .or. &
          size(alpha_y_glob,1) /= nxg_alpha .or. size(alpha_y_glob,2) /= nyg_alpha) then
        cache_local_ok = 0
      endif
    endif
  endif

  ! Require every PE in this RF domain group to agree before reusing root data.
  if (mype == root) then
    cache_all_ok = cache_local_ok
    do n = 2, npelist
      call mpp_recv(peer_cache_ok, glen=1, from_pe=pelist(n), block=.true.)
      cache_all_ok = min(cache_all_ok, peer_cache_ok)
    enddo
    do n = 2, npelist
      call mpp_send(cache_all_ok, plen=1, to_pe=pelist(n))
    enddo
  else
    call mpp_send(cache_local_ok, plen=1, to_pe=root)
    call mpp_recv(cache_all_ok, glen=1, from_pe=root, block=.true.)
  endif
  call mpp_sync_self()

  if (cache_all_ok == 1) then
    deallocate(pelist)
    return
  endif

  ni = iec - isc + 1
  nj = jec - jsc + 1
  n2 = ni * nj
  send_count = 2*n2

  ! -----------------------------
  ! 1) Build local alpha on compute region
  ! -----------------------------
  do j = jsc, jec
    do i = isc, iec
      ax = exp( -dx(i,j) / max(Lx_m, 1.0) )
      ay = exp( -dy(i,j) / max(Ly_m, 1.0) )
      alpha_x(i,j) = min(0.999999, max(0.0, ax))
      alpha_y(i,j) = min(0.999999, max(0.0, ay))
    end do
  end do

  ! -----------------------------
  ! 2) Allocate root global arrays for this active alpha build
  ! -----------------------------
  if (mype == root) then
    if (allocated(alpha_x_glob)) deallocate(alpha_x_glob)
    if (allocated(alpha_y_glob)) deallocate(alpha_y_glob)

    if (skipg == 0) then
      ! Original/full-resolution RF alpha storage.
      allocate(alpha_x_glob(isg:ieg, jsg:jeg))
      allocate(alpha_y_glob(isg:ieg, jsg:jeg))
    else
      ! Store the active coarse alpha only.  The fine alpha arrays are
      !  temporary and are deallocated before leaving this routine.
      allocate(alpha_x_glob(1:nxg_alpha, 1:nyg_alpha))
      allocate(alpha_y_glob(1:nxg_alpha, 1:nyg_alpha))
      allocate(alpha_x_fine_root(isg:ieg, jsg:jeg))
      allocate(alpha_y_fine_root(isg:ieg, jsg:jeg))
      alpha_x_fine_root(:,:) = 0.0
      alpha_y_fine_root(:,:) = 0.0
    endif

    alpha_x_glob(:,:) = 0.0
    alpha_y_glob(:,:) = 0.0
  end if

  ! -----------------------------
  ! 3) Pack local compute slab into sendbuf:
  !    sendbuf(1:n2)       = alpha_x
  !    sendbuf(n2+1:2*n2)  = alpha_y
  ! -----------------------------
  allocate(sendbuf(send_count))

  n = 0
  do j = jsc, jec
    do i = isc, iec
      n = n + 1
      sendbuf(n)       = alpha_x(i,j)
      sendbuf(n + n2)  = alpha_y(i,j)
    end do
  end do

  ! -----------------------------
  ! 4) Root gathers metadata (compute-domain bounds + slab sizes)
  ! -----------------------------
  Ldims(1) = isc
  Ldims(2) = iec
  Ldims(3) = jsc
  Ldims(4) = jec

  if (mype == root) then
    allocate(Gdims(4*npelist))
    allocate(LsizeS(npelist))
    allocate(displ(npelist))

    ! root itself
    Gdims(1:4)  = Ldims(1:4)
    LsizeS(1)   = send_count

    do n = 2, npelist
      call mpp_recv(Gdims((n-1)*4+1), glen=4, from_pe=pelist(n), block=.true.)
      call mpp_recv(LsizeS(n),        glen=1, from_pe=pelist(n), block=.true.)
    end do

    displ(1) = 0
    tot = LsizeS(1)
    do n = 2, npelist
      displ(n) = tot
      tot      = tot + LsizeS(n)
    end do

    allocate(recvbuf(tot))

    ! copy root's own slab first
    recvbuf(1:LsizeS(1)) = sendbuf(1:LsizeS(1))

    ! then receive the rest slabs into the right offsets
    do n = 2, npelist
      call mpp_recv(recvbuf(displ(n)+1), glen=LsizeS(n), from_pe=pelist(n), block=.true.)
    end do

  else
    ! non-root sends meta to root, then slab
    call mpp_send(Ldims(1),   plen=4,    to_pe=root)
    call mpp_send(send_count, plen=1, to_pe=root)
    call mpp_send(sendbuf(1), plen=send_count, to_pe=root)
  end if

  call mpp_sync_self()

  ! -----------------------------
  ! 5) Root reconstructs alpha_x_glob / alpha_y_glob
  ! -----------------------------
  if (mype == root) then
    if (skipg == 0) then
      alpha_x_glob(:,:) = 0.0
      alpha_y_glob(:,:) = 0.0
    else
      alpha_x_fine_root(:,:) = 0.0
      alpha_y_fine_root(:,:) = 0.0
    endif

    do p = 1, npelist
      i1  = Gdims((p-1)*4 + 1)
      i2  = Gdims((p-1)*4 + 2)
      j1  = Gdims((p-1)*4 + 3)
      j2  = Gdims((p-1)*4 + 4)

      nij = (i2 - i1 + 1) * (j2 - j1 + 1)
      off = displ(p)
      base = off

      ! unpack alpha_x first
      kk = 0
      do j = j1, j2
        do i = i1, i2
          kk = kk + 1
          if (skipg == 0) then
            alpha_x_glob(i,j) = recvbuf(base + kk)
          else
            alpha_x_fine_root(i,j) = recvbuf(base + kk)
          endif
        end do
      end do

      ! unpack alpha_y next
      base = off + nij
      kk = 0
      do j = j1, j2
        do i = i1, i2
          kk = kk + 1
          if (skipg == 0) then
            alpha_y_glob(i,j) = recvbuf(base + kk)
          else
            alpha_y_fine_root(i,j) = recvbuf(base + kk)
          endif
        end do
      end do

    end do

    if (skipg > 0) then
      call build_coarse_alpha_from_fine(alpha_x_fine_root, alpha_y_fine_root, &
                                        alpha_x_glob, alpha_y_glob, &
                                        isg, ieg, jsg, jeg, nxg_alpha, nyg_alpha, strideg)
      deallocate(alpha_x_fine_root, alpha_y_fine_root)
    endif

    deallocate(Gdims, LsizeS, displ, recvbuf)
  end if

  deallocate(sendbuf)
  deallocate(pelist)

  rf_isg = isg
  rf_ieg = ieg
  rf_jsg = jsg
  rf_jeg = jeg
  rf_alpha_stride = strideg
  rf_alpha_Lx_m_stored = Lx_m
  rf_alpha_Ly_m_stored = Ly_m
  rf_alpha_mode_stored = 1
  rf_alpha_tile_stored = nudge_tile_num
  rf_alpha_root_stored = root
  rf_alpha_isc_stored = isc
  rf_alpha_iec_stored = iec
  rf_alpha_jsc_stored = jsc
  rf_alpha_jec_stored = jec
  rf_alpha_dx_min_stored = dx_min
  rf_alpha_dx_max_stored = dx_max
  rf_alpha_dy_min_stored = dy_min
  rf_alpha_dy_max_stored = dy_max
  rf_alpha_dx_sum_stored = dx_sum
  rf_alpha_dy_sum_stored = dy_sum
  rf_alpha_dx_sq_sum_stored = dx_sq_sum
  rf_alpha_dy_sq_sum_stored = dy_sq_sum
  rf_alpha_stored = .true.
end subroutine build_rf_alpha_from_L

!> 2-D recursive filter for scalar fields (T, q, etc.)
!!  Conceptual replacement for del2_scalar.
subroutine rf2d_scalar(qdt, kmd, bd, nf, npx, npy, domain, do_norm)
  use mpp_domains_mod, only: domain2d, mpp_update_domains, mpp_get_compute_domain,mpp_get_global_domain
  use mpp_mod,         only: mpp_root_pe
  implicit none

  type(domain2d), intent(inout) :: domain
  integer,        intent(in)    :: kmd, npx, npy, nf
  logical,        intent(in)    :: do_norm

  ! local compute bounds (assumed available in host scope like your file does)
  type(fv_grid_bounds_type), intent(in) :: bd
  real, intent(inout) :: qdt    (is:ie, js:je, kmd)

  real :: norm2d(is:ie, js:je)
  integer :: isg, ieg, jsg, jeg

  ! 1) (Optional) compute norm2d (global RF applied to ones, then gather back)
  if (do_norm) then
     call rf2d_norm_global(norm2d, bd, npx, npy, nf, domain)
  else
     norm2d = 1.0
  end if

  ! 2) Apply global RF to qdt itself
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
  call rf2d_apply_global(qdt, bd, nf, isg, ieg, jsg, jeg, kmd, domain)

  ! 3) Normalize (only meaningful when do_norm=.true.)
  if (do_norm) then
     qdt(is:ie, js:je, 1:kmd) = qdt(is:ie, js:je, 1:kmd) / &
                                           spread(norm2d, dim=3, ncopies=kmd)
  end if

end subroutine rf2d_scalar

!> 2-D recursive filter for vector winds (du, dv).
!!  Conceptual replacement for del2_uv using rf2d_scalar-like smoother.
!------------------------------------------------------------
! Helpers: NaN check (compute-domain only)
!------------------------------------------------------------

subroutine rf_uv(du, dv, kmd, bd, nf, npx, npy, gridstruct, domain)
  use fv_arrays_mod,    only: fv_grid_type, fv_grid_bounds_type, R_GRID
  use mpp_domains_mod,  only: domain2d, mpp_get_compute_domain
  implicit none

  type(domain2d),             intent(inout) :: domain
  type(fv_grid_bounds_type),  intent(in)    :: bd
  type(fv_grid_type),         intent(in), target :: gridstruct
  integer,                    intent(in)    :: kmd, npx, npy, nf

  real, intent(inout) :: du(is:ie, js:je, kmd)
  real, intent(inout) :: dv(is:ie, js:je, kmd)

  real(kind=R_GRID), pointer :: vlon(:,:,:), vlat(:,:,:)
  real, allocatable :: v1(:,:,:), v2(:,:,:), v3(:,:,:)
  real :: norm2d(is:ie,js:je)
  integer :: i,j,k,isg, ieg, jsg, jeg

  vlon => gridstruct%vlon
  vlat => gridstruct%vlat

  allocate(v1(bd%isd:bd%ied, bd%jsd:bd%jed, kmd))
  allocate(v2(bd%isd:bd%ied, bd%jsd:bd%jed, kmd))
  allocate(v3(bd%isd:bd%ied, bd%jsd:bd%jed, kmd))
  v1 = 0.0; v2 = 0.0; v3 = 0.0

  ! --- transform (u,v) -> Cartesian components (v1,v2,v3) on compute domain
!$omp parallel do private(i,j,k) default(shared)
  do k = 1, kmd
    do j = js, je
      do i = is, ie
        v1(i,j,k) = du(i,j,k)*vlon(i,j,1) + dv(i,j,k)*vlat(i,j,1)
        v2(i,j,k) = du(i,j,k)*vlon(i,j,2) + dv(i,j,k)*vlat(i,j,2)
        v3(i,j,k) = du(i,j,k)*vlon(i,j,3) + dv(i,j,k)*vlat(i,j,3)
      end do
    end do
  end do
!$omp end parallel do

  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
  ! --- apply GLOBAL RF.  The variable-coefficient forward/adjoint pair is
  ! symmetric but does not preserve a constant field without normalization.
  ! Compute that normalization once and reuse it for all three Cartesian
  ! wind components.
  call rf2d_norm_global(norm2d, bd, npx, npy, nf, domain)
  call rf2d_apply_global(v1(is:ie, js:je, 1:kmd), bd, nf, isg, ieg, jsg, jeg, kmd, domain)
  call rf2d_apply_global(v2(is:ie, js:je, 1:kmd), bd, nf, isg, ieg, jsg, jeg, kmd, domain)
  call rf2d_apply_global(v3(is:ie, js:je, 1:kmd), bd, nf, isg, ieg, jsg, jeg, kmd, domain)
  v1(is:ie,js:je,1:kmd)=v1(is:ie,js:je,1:kmd)/spread(norm2d,dim=3,ncopies=kmd)
  v2(is:ie,js:je,1:kmd)=v2(is:ie,js:je,1:kmd)/spread(norm2d,dim=3,ncopies=kmd)
  v3(is:ie,js:je,1:kmd)=v3(is:ie,js:je,1:kmd)/spread(norm2d,dim=3,ncopies=kmd)

  ! --- inverse transform back to (u,v)
!$omp parallel do private(i,j,k) default(shared)
  do k = 1, kmd
    do j = js, je
      do i = is, ie
        du(i,j,k) = v1(i,j,k)*vlon(i,j,1) + v2(i,j,k)*vlon(i,j,2) + v3(i,j,k)*vlon(i,j,3)
        dv(i,j,k) = v1(i,j,k)*vlat(i,j,1) + v2(i,j,k)*vlat(i,j,2) + v3(i,j,k)*vlat(i,j,3)
      end do
    end do
  end do
!$omp end parallel do

  deallocate(v1, v2, v3)
end subroutine rf_uv

!=========================================================
! Improved rf2d_apply_global
!=========================================================
subroutine rf2d_apply_global(field_inout, bd, nf, isg, ieg, jsg, jeg, kmd, domain)
  use mpp_mod,         only: mpp_pe, mpp_root_pe, mpp_send, mpp_recv, mpp_sync_self, &
                              mpp_npes, mpp_broadcast, mpp_error, FATAL
  use mpp_domains_mod, only: domain2d, mpp_get_pelist, mpp_update_domains
  use fv_arrays_mod,   only: fv_grid_bounds_type
  implicit none

  type(domain2d),            intent(inout) :: domain
  type(fv_grid_bounds_type), intent(in)    :: bd
  integer,                   intent(in)    :: nf, isg, ieg, jsg, jeg, kmd
  real,                      intent(inout) :: field_inout(is:ie, js:je, 1:kmd)

  integer :: mype, root, npelist, nact_loc
  integer :: ni, nj, ni_skipg, nj_skipg, local_send_count
  integer :: i, j, k, n, idx
  integer :: skipg, strideg
  integer :: i_global, j_global, ii, jj, ii_next, jj_next
  integer :: i_coarse_global, j_coarse_global
  real    :: wi, wj
  logical :: is_root
  real    :: field(bd%isd:bd%ied, bd%jsd:bd%jed, 1:kmd)

  integer, allocatable :: pelist(:)
  integer, allocatable :: pe_is(:), pe_ie(:), pe_js(:), pe_je(:)
  integer, allocatable :: pe_sizes(:)

  real, allocatable :: send_buf(:)
  real, allocatable :: recv_buf(:)

  real, allocatable :: global_field(:,:,:)
  real, allocatable :: global_field_coarse(:,:,:)
  ! Coarse alpha is now built/stored directly in build_rf_alpha_from_L
  !  when nf_skip_grid > 0, so no temporary coarse-alpha arrays are needed here.

  integer :: nxg_coarse, nyg_coarse
  integer :: is_coarse, ie_coarse, js_coarse, je_coarse
  integer :: first_i_skipg, last_i_skipg, first_j_skipg, last_j_skipg
  integer :: i_last_skipg, j_last_skipg

  ! Comparison switch for option 4.  When disabled, fall through to the
  !  original per-tile gather/filter/scatter implementation below.
  if (nudge_option == 4 .and. rf_cube_continuous) then
    call rf2d_apply_cube_global(field_inout, bd, nf, isg, ieg, jsg, jeg, kmd, domain)
    return
  endif

  mype = mpp_pe()
  call get_rf_domain_pelist(domain, pelist, npelist, root)
  is_root = (mype == root)

  skipg = max(0, nf_skip_grid)
  strideg = skipg + 1

  nact_loc = min(nact, kmd)

  !$OMP parallel do default(none) shared(is,ie,js,je,kmd,field,field_inout)
  do k = 1, kmd
    do j = js, je
      do i = is, ie
        field(i,j,k) = field_inout(i,j,k)
      enddo
    enddo
  enddo
  !$OMP end parallel do

  call mpp_update_domains(field, domain, complete=.true.)

  ni = ie - is + 1
  nj = je - js + 1

  ! =================================================================
  ! CASE 1: NO STRIDE
  ! =================================================================
  if (skipg == 0) then

    if (is_root) then
      allocate(global_field(isg:ieg, jsg:jeg, 1:nact_loc))
      global_field = 0.0
    end if

    if (is_root) then
      allocate(pe_is(npelist), pe_ie(npelist), pe_js(npelist), pe_je(npelist))
      allocate(pe_sizes(npelist))

      pe_is(1) = is
      pe_ie(1) = ie
      pe_js(1) = js
      pe_je(1) = je
      local_send_count = ni * nj * nact_loc
      pe_sizes(1) = local_send_count

      do n = 2, npelist
        call mpp_recv(pe_is(n),    glen=1, from_pe=pelist(n), block=.true.)
        call mpp_recv(pe_ie(n),    glen=1, from_pe=pelist(n), block=.true.)
        call mpp_recv(pe_js(n),    glen=1, from_pe=pelist(n), block=.true.)
        call mpp_recv(pe_je(n),    glen=1, from_pe=pelist(n), block=.true.)
        call mpp_recv(pe_sizes(n), glen=1, from_pe=pelist(n), block=.true.)
      enddo
    else
      call mpp_send(is, plen=1, to_pe=root)
      call mpp_send(ie, plen=1, to_pe=root)
      call mpp_send(js, plen=1, to_pe=root)
      call mpp_send(je, plen=1, to_pe=root)
      local_send_count = ni * nj * nact_loc
      call mpp_send(local_send_count, plen=1, to_pe=root)
    endif
    call mpp_sync_self()

    local_send_count = ni * nj * nact_loc
    allocate(send_buf(local_send_count))

    idx = 1
    do k = 1, nact_loc
      do j = js, je
        do i = is, ie
          send_buf(idx) = field(i, j, klist(k))
          idx = idx + 1
        enddo
      enddo
    enddo

    if (is_root) then
      allocate(recv_buf(maxval(pe_sizes)))

      do n = 1, npelist
        if (pelist(n) == root) then
          recv_buf(1:pe_sizes(n)) = send_buf(1:pe_sizes(n))
        else
          call mpp_recv(recv_buf(1), glen=pe_sizes(n), from_pe=pelist(n), block=.true.)
        endif

        idx = 1
        do k = 1, nact_loc
          do j = pe_js(n), pe_je(n)
            do i = pe_is(n), pe_ie(n)
              global_field(i,j,k) = recv_buf(idx)
              idx = idx + 1
            enddo
          enddo
        enddo
      enddo

      deallocate(recv_buf)
    else
      call mpp_send(send_buf(1), plen=local_send_count, to_pe=root)
    endif
    call mpp_sync_self()

    deallocate(send_buf)

    if (is_root) then
      !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(STATIC)
      do k = 1, nact_loc
        call rf1d_x_global_2d(global_field(:,:,k), alpha_x_glob, nf, isg, ieg, jsg, jeg)
        call rf1d_y_global_2d(global_field(:,:,k), alpha_y_glob, nf, isg, ieg, jsg, jeg)
      enddo
      !$OMP END PARALLEL DO
    endif

    if (.not. is_root) then
      allocate(global_field(isg:ieg, jsg:jeg, 1:nact_loc))
    endif

    call rf_pelist_broadcast_real(global_field(isg,jsg,1), size(global_field), root, pelist, npelist)

    !$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(is,ie,js,je,nact_loc,field_inout,global_field,klist) &
    !$OMP SCHEDULE(STATIC)
    do k = 1, nact_loc
      do j = js, je
        do i = is, ie
          field_inout(i,j,klist(k)) = global_field(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(global_field)
    if (is_root) deallocate(pe_is, pe_ie, pe_js, pe_je, pe_sizes)
    deallocate(pelist)
    return
  endif

  ! =================================================================
  ! CASE 2: WITH STRIDE
  ! =================================================================
  ni_skipg = 0
  do i = is, ie
    if (mod(i - isg, strideg) == 0) ni_skipg = ni_skipg + 1
  enddo

  nj_skipg = 0
  do j = js, je
    if (mod(j - jsg, strideg) == 0) nj_skipg = nj_skipg + 1
  enddo

  nxg_coarse = (ieg - isg) / strideg + 1
  nyg_coarse = (jeg - jsg) / strideg + 1

  if (mype == root) then
    print *, 'RF with skipg=', skipg
    print *, '  Global domain:', isg, 'to', ieg, 'x', jsg, 'to', jeg
    print *, '  Coarse grid:', nxg_coarse, 'x', nyg_coarse, 'points'
    print *, '  Reduction:', real((ieg-isg+1)*(jeg-jsg+1))/real(nxg_coarse*nyg_coarse), 'x'
  endif

  if (is_root) then
    allocate(global_field_coarse(1:nxg_coarse, 1:nyg_coarse, 1:nact_loc))
    global_field_coarse = 0.0
  endif

  if (is_root) then
    allocate(pe_is(npelist), pe_ie(npelist), pe_js(npelist), pe_je(npelist))
    allocate(pe_sizes(npelist))

    first_i_skipg = 0
    last_i_skipg = 0
    do i = is, ie
      if (mod(i - isg, strideg) == 0) then
        if (first_i_skipg == 0) first_i_skipg = i
        last_i_skipg = i
      endif
    enddo

    first_j_skipg = 0
    last_j_skipg = 0
    do j = js, je
      if (mod(j - jsg, strideg) == 0) then
        if (first_j_skipg == 0) first_j_skipg = j
        last_j_skipg = j
      endif
    enddo

    if (first_i_skipg > 0) then
      pe_is(1) = (first_i_skipg - isg) / strideg + 1
      pe_ie(1) = (last_i_skipg  - isg) / strideg + 1
    else
      pe_is(1) = 1
      pe_ie(1) = 0
    endif

    if (first_j_skipg > 0) then
      pe_js(1) = (first_j_skipg - jsg) / strideg + 1
      pe_je(1) = (last_j_skipg  - jsg) / strideg + 1
    else
      pe_js(1) = 1
      pe_je(1) = 0
    endif

    local_send_count = ni_skipg * nj_skipg * nact_loc
    pe_sizes(1) = local_send_count

    do n = 2, npelist
      call mpp_recv(pe_is(n),    glen=1, from_pe=pelist(n), block=.true.)
      call mpp_recv(pe_ie(n),    glen=1, from_pe=pelist(n), block=.true.)
      call mpp_recv(pe_js(n),    glen=1, from_pe=pelist(n), block=.true.)
      call mpp_recv(pe_je(n),    glen=1, from_pe=pelist(n), block=.true.)
      call mpp_recv(pe_sizes(n), glen=1, from_pe=pelist(n), block=.true.)
    enddo
  else
    first_i_skipg = 0
    last_i_skipg = 0
    do i = is, ie
      if (mod(i - isg, strideg) == 0) then
        if (first_i_skipg == 0) first_i_skipg = i
        last_i_skipg = i
      endif
    enddo

    first_j_skipg = 0
    last_j_skipg = 0
    do j = js, je
      if (mod(j - jsg, strideg) == 0) then
        if (first_j_skipg == 0) first_j_skipg = j
        last_j_skipg = j
      endif
    enddo

    if (first_i_skipg > 0) then
      is_coarse = (first_i_skipg - isg) / strideg + 1
      ie_coarse = (last_i_skipg  - isg) / strideg + 1
    else
      is_coarse = 1
      ie_coarse = 0
    endif

    if (first_j_skipg > 0) then
      js_coarse = (first_j_skipg - jsg) / strideg + 1
      je_coarse = (last_j_skipg  - jsg) / strideg + 1
    else
      js_coarse = 1
      je_coarse = 0
    endif

    call mpp_send(is_coarse, plen=1, to_pe=root)
    call mpp_send(ie_coarse, plen=1, to_pe=root)
    call mpp_send(js_coarse, plen=1, to_pe=root)
    call mpp_send(je_coarse, plen=1, to_pe=root)
    local_send_count = ni_skipg * nj_skipg * nact_loc
    call mpp_send(local_send_count, plen=1, to_pe=root)
  endif
  call mpp_sync_self()

  local_send_count = ni_skipg * nj_skipg * nact_loc
  allocate(send_buf(max(1, local_send_count)))

  if (ni_skipg > 0 .and. nj_skipg > 0) then
    idx = 1
    do k = 1, nact_loc
      do j = js, je
        if (mod(j - jsg, strideg) == 0) then
          do i = is, ie
            if (mod(i - isg, strideg) == 0) then
              send_buf(idx) = field(i, j, klist(k))
              idx = idx + 1
            endif
          enddo
        endif
      enddo
    enddo
  endif

  if (is_root) then
    allocate(recv_buf(max(1, maxval(pe_sizes))))

    do n = 1, npelist
      if (pe_sizes(n) > 0) then
        if (pelist(n) == root) then
          recv_buf(1:pe_sizes(n)) = send_buf(1:pe_sizes(n))
        else
          call mpp_recv(recv_buf(1), glen=pe_sizes(n), from_pe=pelist(n), block=.true.)
        endif

        idx = 1
        do k = 1, nact_loc
          do j = pe_js(n), pe_je(n)
            do i = pe_is(n), pe_ie(n)
              global_field_coarse(i,j,k) = recv_buf(idx)
              idx = idx + 1
            enddo
          enddo
        enddo
      endif
    enddo

    deallocate(recv_buf)
  else
    if (ni_skipg * nj_skipg * nact_loc > 0) then
      call mpp_send(send_buf(1), plen=local_send_count, to_pe=root)
    endif
  endif
  call mpp_sync_self()

  deallocate(send_buf)

  if (is_root) then
    !$OMP PARALLEL DO PRIVATE(k) SCHEDULE(STATIC)
    do k = 1, nact_loc
      call rf1d_x_global_2d(global_field_coarse(:,:,k), alpha_x_glob, nf, &
                            1, nxg_coarse, 1, nyg_coarse)
      call rf1d_y_global_2d(global_field_coarse(:,:,k), alpha_y_glob, nf, &
                            1, nxg_coarse, 1, nyg_coarse)
    enddo
    !$OMP END PARALLEL DO
  endif

  if (.not. is_root) then
    allocate(global_field_coarse(1:nxg_coarse, 1:nyg_coarse, 1:nact_loc))
  endif

  call rf_pelist_broadcast_real(global_field_coarse(1,1,1), size(global_field_coarse), root, pelist, npelist)

  i_last_skipg = isg + (nxg_coarse - 1) * strideg
  j_last_skipg = jsg + (nyg_coarse - 1) * strideg

  !$OMP PARALLEL DO PRIVATE(k,j,i,j_global,i_global,jj,ii,jj_next,ii_next, &
  !$OMP j_coarse_global,i_coarse_global,wj,wi) &
  !$OMP SHARED(js,je,is,ie,nact_loc,field_inout,global_field_coarse, &
  !$OMP        jsg,jeg,isg,ieg,strideg,nxg_coarse,nyg_coarse, &
  !$OMP        i_last_skipg,j_last_skipg,klist) &
  !$OMP SCHEDULE(STATIC)
  do k = 1, nact_loc
    do j = js, je
      j_global = j

      if (j_global > j_last_skipg) then
        jj = nyg_coarse - 1
        jj_next = nyg_coarse
        wj = real(j_global - j_last_skipg) / real(strideg)
      else
        jj = (j_global - jsg) / strideg + 1
        jj = max(1, min(jj, nyg_coarse - 1))
        jj_next = min(jj + 1, nyg_coarse)
        j_coarse_global = jsg + (jj - 1) * strideg
        wj = real(j_global - j_coarse_global) / real(strideg)
        wj = max(0.0, min(1.0, wj))
      endif

      do i = is, ie
        i_global = i

        if (i_global > i_last_skipg) then
          ii = nxg_coarse - 1
          ii_next = nxg_coarse
          wi = real(i_global - i_last_skipg) / real(strideg)
        else
          ii = (i_global - isg) / strideg + 1
          ii = max(1, min(ii, nxg_coarse - 1))
          ii_next = min(ii + 1, nxg_coarse)
          i_coarse_global = isg + (ii - 1) * strideg
          wi = real(i_global - i_coarse_global) / real(strideg)
          wi = max(0.0, min(1.0, wi))
        endif

        field_inout(i,j,klist(k)) = &
          (1.0 - wi) * (1.0 - wj) * global_field_coarse(ii,      jj,      k) + &
          wi         * (1.0 - wj) * global_field_coarse(ii_next, jj,      k) + &
          (1.0 - wi) * wj         * global_field_coarse(ii,      jj_next, k) + &
          wi         * wj         * global_field_coarse(ii_next, jj_next, k)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  deallocate(global_field_coarse)
  if (is_root) deallocate(pe_is, pe_ie, pe_js, pe_je, pe_sizes)
  deallocate(pelist)

end subroutine rf2d_apply_global

subroutine rf1d_x_global_2d(f, alpha2d, nf, isg, ieg, jsg, jeg)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  integer, intent(in) :: nf, isg, ieg, jsg, jeg
  real,    intent(inout) :: f(isg:ieg, jsg:jeg)
  real,    intent(in)    :: alpha2d(isg:ieg, jsg:jeg)

  integer :: it, i, j
  real :: a, b
  real, allocatable :: tmp(:)

  allocate(tmp(isg:ieg))

  do it = 1, nf
    do j = jsg, jeg
      tmp(isg) = f(isg,j)
      do i = isg+1, ieg
        a = alpha2d(i,j)
        if (ieee_is_nan(a)) a = 0.0
        if (a < 0.0) a = 0.0
        if (a > 0.999999) a = 0.999999
        b = 1.0 - a
        tmp(i) = a*tmp(i-1) + b*f(i,j)
      enddo
      ! Exact adjoint of the forward one-sided recursion.  First solve
      ! G^T without the local (1-alpha) factor, then apply that factor at
      ! each point that used it in the forward pass.  The left endpoint has
      ! forward coefficient 1 and therefore is not multiplied by (1-alpha).
      f(ieg,j) = tmp(ieg)
      do i = ieg-1, isg, -1
        a = alpha2d(i+1,j)
        if (ieee_is_nan(a)) a = 0.0
        if (a < 0.0) a = 0.0
        if (a > 0.999999) a = 0.999999
        f(i,j) = tmp(i) + a*f(i+1,j)
      enddo
      do i = isg+1, ieg
        a = alpha2d(i,j)
        if (ieee_is_nan(a)) a = 0.0
        if (a < 0.0) a = 0.0
        if (a > 0.999999) a = 0.999999
        f(i,j) = (1.0-a)*f(i,j)
      enddo
    enddo
  enddo

  deallocate(tmp)
end subroutine rf1d_x_global_2d

subroutine rf1d_y_global_2d(f, alpha2d, nf, isg, ieg, jsg, jeg)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  integer, intent(in) :: nf, isg, ieg, jsg, jeg
  real,    intent(inout) :: f(isg:ieg, jsg:jeg)
  real,    intent(in)    :: alpha2d(isg:ieg, jsg:jeg)

  integer :: it, i, j
  real :: a, b
  real, allocatable :: tmp(:)

  allocate(tmp(jsg:jeg))

  do it = 1, nf
    do i = isg, ieg
      tmp(jsg) = f(i,jsg)
      do j = jsg+1, jeg
        a = alpha2d(i,j)
        if (ieee_is_nan(a)) a = 0.0
        if (a < 0.0) a = 0.0
        if (a > 0.999999) a = 0.999999
        b = 1.0 - a
        tmp(j) = a*tmp(j-1) + b*f(i,j)
      enddo
      ! Exact adjoint of the forward one-sided recursion; see the X pass.
      f(i,jeg) = tmp(jeg)
      do j = jeg-1, jsg, -1
        a = alpha2d(i,j+1)
        if (ieee_is_nan(a)) a = 0.0
        if (a < 0.0) a = 0.0
        if (a > 0.999999) a = 0.999999
        f(i,j) = tmp(j) + a*f(i,j+1)
      enddo
      do j = jsg+1, jeg
        a = alpha2d(i,j)
        if (ieee_is_nan(a)) a = 0.0
        if (a < 0.0) a = 0.0
        if (a > 0.999999) a = 0.999999
        f(i,j) = (1.0-a)*f(i,j)
      enddo
    enddo
  enddo

  deallocate(tmp)
end subroutine rf1d_y_global_2d

subroutine rf2d_norm_global(norm2d, bd, npx, npy, nf, domain)
  use mpp_domains_mod, only: domain2d,mpp_get_global_domain
  use fv_arrays_mod,   only: fv_grid_bounds_type
  implicit none

  type(domain2d),            intent(inout) :: domain
  type(fv_grid_bounds_type), intent(in)    :: bd
  integer,                  intent(in)     :: nf, npx, npy
  real,                     intent(out)    :: norm2d(is:ie, js:je)

  real, allocatable :: ones(:,:,:)
  integer :: nact_save
  integer, allocatable :: klist_save(:)
  logical :: klist_was_alloc
  integer :: kmd,isg,ieg,jsg,jeg

  kmd = 1
  allocate(ones(is:ie, js:je, 1:kmd))
  ones(:,:,1) = 1.0

  ! ---- Save current module state (nact/klist) ----
  nact_save = nact
  klist_was_alloc = allocated(klist)
  if (klist_was_alloc) then
     allocate(klist_save(size(klist)))
     klist_save(:) = klist(:)
  else
     allocate(klist_save(0))
  end if

  ! ---- Override for norm call: only 1 active level ----
  nact = 1
  if (allocated(klist)) deallocate(klist)
  allocate(klist(1))
  klist(1) = 1

  ! ---- Call RF on ones (3D array, kmd=1) ----
  call mpp_get_global_domain(domain, isg, ieg, jsg, jeg)
  call rf2d_apply_global(ones, bd, nf, isg, ieg, jsg, jeg, 1, domain)

  norm2d(:,:) = ones(:,:,1)

  ! ---- Restore module state ----
  nact = nact_save
  if (allocated(klist)) deallocate(klist)
  if (klist_was_alloc) then
     allocate(klist(size(klist_save)))
     klist(:) = klist_save(:)
  end if
  if (allocated(klist_save)) deallocate(klist_save)

  ! safer than array-temp MAX if you saw "stack smashing" before:
  where (norm2d(:,:) < 1.0e-12) norm2d(:,:) = 1.0e-12

  deallocate(ones)
end subroutine rf2d_norm_global

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Adapted from GSI!!!!!!!!!!!!!!!!!!!!!
 subroutine rmse_bias(a, rms, bias, area)
   real, intent(in):: a(is:ie,js:je)
   real, intent(IN) :: area(isd:ied,jsd:jed)
   real, intent(out):: rms, bias
   integer:: i,j
   real:: total_area

   total_area = 4.*pi*radius**2

    rms = 0.
   bias = 0.
   do j=js,je
      do i=is,ie
         bias = bias + area(i,j) * a(i,j)
          rms = rms  + area(i,j) * a(i,j)**2
      enddo
   enddo
   call mp_reduce_sum(bias)
   call mp_reduce_sum(rms)

   bias = bias / total_area
    rms = sqrt( rms / total_area )
 end subroutine rmse_bias

 subroutine corr(a, b, co, area)
 real, intent(in):: a(is:ie,js:je), b(is:ie,js:je)
 real, intent(in):: area(isd:ied,jsd:jed)
 real, intent(out):: co
 real:: m_a, m_b, std_a, std_b
 integer:: i,j
 real:: total_area

   total_area = 4.*pi*radius**2

! Compute standard deviation:
   call std(a, m_a, std_a, area)
   call std(b, m_b, std_b, area)

! Compute correlation:
   co = 0.
   do j=js,je
      do i=is,ie
         co = co + area(i,j) * (a(i,j)-m_a)*(b(i,j)-m_b)
      enddo
   enddo
   call mp_reduce_sum(co)
   co = co / (total_area*std_a*std_b )
 end subroutine corr

 subroutine std(a, mean, stdv, area)
 real,  intent(in):: a(is:ie,js:je)
 real, intent(IN) :: area(isd:ied,jsd:jed)
 real, intent(out):: mean, stdv
 integer:: i,j
 real:: total_area

   total_area = 4.*pi*radius**2

   mean = 0.
   do j=js,je
      do i=is,ie
         mean = mean + area(i,j) * a(i,j)
      enddo
   enddo
   call mp_reduce_sum(mean)
   mean = mean / total_area

   stdv = 0.
   do j=js,je
      do i=is,ie
         stdv = stdv + area(i,j) * (a(i,j)-mean)**2
      enddo
   enddo
   call mp_reduce_sum(stdv)
   stdv = sqrt( stdv / total_area )
 end subroutine std

end module fv_nwp_nudge_mod
