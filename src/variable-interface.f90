!    A suite of routines to interpolate Model for Prediction Across
!    Scales (MPAS) Atmosphere (MPAS-A) core variables to a Gaussian
!    grid (of user specified spectral truncation) and readable by the
!    Gridded Analysis and Display System (GrADS).

!    Copyright (C) 2013 Henry R. Winterbottom

!    Email: Henry.Winterbottom@noaa.gov

!    Snail-mail:

!    Henry R. Winterbottom
!    NOAA/OAR/PSD R/PSD1
!    325 Broadway
!    Boulder, CO 80303-3328

!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation; either version 2 of
!    the License, or (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!    General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!    02110-1301 USA.

module variable_interface

  !=======================================================================
  
  ! Define associated modules and subroutines
  
  !-----------------------------------------------------------------------
  
  use constants
  use kinds

  !-----------------------------------------------------------------------

  use interpolation_interface
  use mpi_interface
  use namelist
  use netcdf

  !-----------------------------------------------------------------------
  
  implicit none

  !-----------------------------------------------------------------------

  ! Define all data and structure types for routine; these variables
  ! are variables required by the subroutines within this module

  type variable_info
     character(len=50)                                        :: grads_string
     character(len=30)                                        :: mpas_name
     character(len=20)                                        :: mpas_hori_coord
     character(len=20)                                        :: mpas_vert_coord
     character(len=10)                                        :: grads_id
     character(len=2)                                         :: interp_type
     logical                                                  :: is_profile
     logical                                                  :: is_derived
     real(r_kind),               dimension(:),    allocatable :: array
     integer                                                  :: xdim
     integer                                                  :: ydim
     integer                                                  :: zdim
     integer                                                  :: tdim
  end type variable_info

  ! Define global variables

  real(r_kind)                                                :: spval             = 1.e30
  integer                                                     :: mpas_ncoords
  integer                                                     :: mpas_nlevels
  integer                                                     :: var_ncoords
  integer                                                     :: var_nlevels
  integer                                                     :: number_of_var
  integer                                                     :: var_counter
  integer                                                     :: ncfileid
  integer                                                     :: ncvarid
  integer                                                     :: ncdimid
  integer                                                     :: ncstatus

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! variable_interface_num_var.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_num_var()

    !=====================================================================

    ! Initialize local variable

    number_of_var = 0

    !---------------------------------------------------------------------

    ! Update local variable accordingly

    if(is_psfc)   number_of_var = number_of_var + 1
    if(is_pslp)   number_of_var = number_of_var + 1
    if(is_u10m)   number_of_var = number_of_var + 1
    if(is_v10m)   number_of_var = number_of_var + 1
    if(is_t2m)    number_of_var = number_of_var + 1
    if(is_sst)    number_of_var = number_of_var + 1
    if(is_ter)    number_of_var = number_of_var + 1
    if(is_uwnd)   number_of_var = number_of_var + 1
    if(is_vwnd)   number_of_var = number_of_var + 1
    if(is_vort)   number_of_var = number_of_var + 1
    if(is_divg)   number_of_var = number_of_var + 1
    if(is_rh)     number_of_var = number_of_var + 1
    if(is_qv)     number_of_var = number_of_var + 1
    if(is_temp)   number_of_var = number_of_var + 1
    if(is_vtemp)  number_of_var = number_of_var + 1
    if(is_theta)  number_of_var = number_of_var + 1
    if(is_thetae) number_of_var = number_of_var + 1
    if(is_psi)    number_of_var = number_of_var + 1
    if(is_chi)    number_of_var = number_of_var + 1
    if(is_pwat)   number_of_var = number_of_var + 1
    if(is_pv)     number_of_var = number_of_var + 1
    if(is_cori)   number_of_var = number_of_var + 1
    if(is_dbz)    number_of_var = number_of_var + 1
    if(is_cdbz)   number_of_var = number_of_var + 1
    if(is_ght)    number_of_var = number_of_var + 1

    !=====================================================================

  end subroutine variable_interface_num_var

  !=======================================================================

  ! variable_interface_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_initialize(grid)

    ! Define variables computed within routine

    type(variable_info),     dimension(:),       allocatable :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid)) allocate(grid(number_of_var))

    !=====================================================================

  end subroutine variable_interface_initialize

  !=======================================================================

  ! variable_interface_cleanup.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_cleanup(grid)

    ! Define variables computed within routine

    type(variable_info),     dimension(:),       allocatable :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(allocated(grid)) deallocate(grid)

    !=====================================================================

  end subroutine variable_interface_cleanup

  !=======================================================================

  ! variable_interface_isolevel.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_isolevel(grid)

    ! Define variable returned by routine

    type(variable_info),     dimension(number_of_var)        :: grid

    !=====================================================================

    ! Define local variables accordingly

    if(is_psfc) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'surface pressure (Pa)'
       grid(var_counter)%mpas_name       = 'surface_pressure'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'psfc'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = 1

    end if ! if(is_psfc)

    ! Define local variables accordingly

    if(is_u10m) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'zonal 10-m wind (m/s)'
       grid(var_counter)%mpas_name       = 'u10'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'u10m'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = 1

    end if ! if(is_u10m)

    ! Define local variables accordingly

    if(is_v10m) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'meridional 10-m wind (m/s)'
       grid(var_counter)%mpas_name       = 'v10'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'v10m'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = 1

    end if ! if(is_v10m)

    ! Define local variables accordingly

    if(is_t2m) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = '2-meter temperature (K)'
       grid(var_counter)%mpas_name       = 't2m'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 't2m'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = 1

    end if ! if(is_t2m)

    ! Define local variables accordingly

    if(is_sst) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'sea-surface/skin temperature (K)'
       grid(var_counter)%mpas_name       = 'sst'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'sst'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = 1

    end if ! if(is_sst)

    ! Define local variables accordingly

    if(is_pwat) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'precipitable water (cm)'
       grid(var_counter)%mpas_name       = 'precipw'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'pwat'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = 1

    end if ! if(is_pwat)

    !=====================================================================

  end subroutine variable_interface_isolevel

  !=======================================================================

  ! variable_interface_profile.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_profile(grid)

    ! Define variable returned by routine

    type(variable_info),     dimension(number_of_var)        :: grid

    !=====================================================================

    ! Define local variable accordingly

    if((.not. is_pinterp) .and. (.not. is_zinterp) .and. (.not.            &
         & is_teinterp)) then

       ! Open external file

       ncstatus = nf90_open(path=trim(variable_filename),                  &
            & mode=nf90_nowrite,ncid=ncfileid)
    
       ! Define local variables

       ncstatus = nf90_inq_dimid(ncfileid,'nVertLevels',ncdimid)
       ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,len=nlev)

       ! Close external file

       ncstatus = nf90_close(ncfileid)

    end if ! if((.not. is_pinterp) .and. (.not. is_zinterp) .and. (.not.   &
           ! is_teinterp))

    !---------------------------------------------------------------------

    ! Define local variables accordingly

    if(is_uwnd) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'zonal wind (m/s)'
       grid(var_counter)%mpas_name       = 'uReconstructZonal'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'u'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_uwnd)

    ! Define local variables accordingly

    if(is_vwnd) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'meridional wind (m/s)'
       grid(var_counter)%mpas_name       = 'uReconstructMeridional'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'v'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_vwnd)

    ! Define local variables accordingly

    if(is_theta) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'potential temperature (K)'
       grid(var_counter)%mpas_name       = 'theta'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'theta'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_theta)

    ! Define local variables accordingly

    if(is_qv) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'water vapor mixing ratio (kg/kg)'
       grid(var_counter)%mpas_name       = 'qv'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'qv'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .false.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_qv)

    !=====================================================================

  end subroutine variable_interface_profile

  !=======================================================================

  ! variable_interface_derived_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_derived_initialize(grid)

    ! Define variable returned by routine

    type(variable_info),     dimension(number_of_var)        :: grid

    !=====================================================================

    ! Define local variables accordingly

    if(is_ter) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'terrain elevation (m)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'terrain'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = 1

    end if ! if(is_ter)

    ! Define local variables accordingly

    if(is_pslp) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'sea-level pressure (Pa)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'pslp'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = 1

    end if ! if(is_pslp)

    ! Define local variables accordingly

    if(is_cori) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'Coriolis parameter (1/s)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'cori'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = 1

    end if ! if(is_cori)

    ! Define local variables accordingly

    if(is_cdbz) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'composite reflectivity (dbz)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%grads_id        = 'cdbz'
       grid(var_counter)%is_profile      = .false.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = 1

    end if ! if(is_cdbz)

    ! Define local variables accordingly

    if(is_vort) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'relative vorticity (1/s)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'vort'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_vort)

    ! Define local variables accordingly

    if(is_divg) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'divergence (1/s)'
       grid(var_counter)%mpas_name       = 'divergence'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'divg'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_divg)

    ! Define local variables accordingly

    if(is_rh) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'relative humidity (%)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'rh'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_rh)

    ! Define local variables accordingly

    if(is_temp) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'temperature (K)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 't'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_temp)

    ! Define local variables accordingly

    if(is_vtemp) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'virtual temperature (K)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'tv'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_vtemp)

    ! Define local variables accordingly

    if(is_thetae) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'equivalent potential temperature (K)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'thetae'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_thetae)

    ! Define local variables accordingly

    if(is_psi) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'streamfunction (m^2/s)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'psi'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_psi)

    ! Define local variables accordingly

    if(is_chi) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'velocity potential (m^2/s)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'chi'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_chi)

    ! Define local variables accordingly

    if(is_pv) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'potential vorticity (PVU)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'pv'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_pv)

    ! Define local variables accordingly

    if(is_dbz) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'simulated radar reflectivity (dbz)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'dbz'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_dbz)

    ! Define local variables accordingly

    if(is_ght) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%interp_type     = 'ba'
       grid(var_counter)%grads_string    = 'geopotential (gpm)'
       grid(var_counter)%mpas_hori_coord = 'nCells'
       grid(var_counter)%mpas_vert_coord = 'nVertLevels'
       grid(var_counter)%grads_id        = 'ght'
       grid(var_counter)%is_profile      = .true.
       grid(var_counter)%is_derived      = .true.
       grid(var_counter)%zdim            = nlev

    end if ! if(is_ght)

    !=====================================================================

  end subroutine variable_interface_derived_initialize

  !=======================================================================

  ! variable_interface_process.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_process(variable_info_grid,srcgrid,        &
       & dstgrid,grid)

    ! Define variable returned by routine

    type(variable_info),     dimension(number_of_var)              :: variable_info_grid
    type(interpgrid)                                               :: srcgrid
    type(interpgrid)                                               :: dstgrid
    type(interppres)                                               :: srcgrid_pres
    type(interppres)                                               :: dstgrid_pres
    type(interphght)                                               :: srcgrid_hght
    type(interphght)                                               :: dstgrid_hght
    type(interpthetae)                                             :: srcgrid_thee
    type(interpthetae)                                             :: dstgrid_thee
    type(grid_interface)                                           :: grid

    ! Define variables computed within routine

    real(r_double),          dimension(:,:),           allocatable :: workgrid_2d
    real(r_double),          dimension(:),             allocatable :: workgrid_1d
    real(r_kind),            dimension(:,:),           allocatable :: workgrid_2r
    real(r_kind),            dimension(:),             allocatable :: workgrid_1r
    integer                                                        :: hori_dim
    integer                                                        :: vert_dim
    integer                                                        :: array_start
    integer                                                        :: array_end

    ! Define counting variables

    integer                                                        :: i, j, k, l
    integer                                                        :: count

    !=====================================================================

    ! Define local variables

    srcgrid%ncoords   = grid%atmos_mass_ncells
    srcgrid%nvertlevs = mpas_nlevels
    srcgrid%npasses   = barnes_npasses
    srcgrid%neighbors = barnes_nneighbors
    dstgrid%ncoords   = grid%gauss_mass_ncells
    dstgrid%nvertlevs = mpas_nlevels
    dstgrid%npasses   = barnes_npasses
    dstgrid%neighbors = barnes_nneighbors

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(srcgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%nvertlevs,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%npasses,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%nvertlevs,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%npasses,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)

    ! Initialize local variables

    call interpolation_initialize_grid(srcgrid)
    call interpolation_initialize_grid(dstgrid)
    call interpolation_initialize_task_balance(srcgrid)
    call interpolation_initialize_task_balance(dstgrid)

    ! Define local variables

    srcgrid%xlong    = real(grid%atmos_mass_lon)
    srcgrid%xlat     = real(grid%atmos_mass_lat)
    srcgrid%weights  = grid%gauss_mass_2_atmos_mass_weights
    srcgrid%grdnbors = grid%gauss_mass_2_atmos_mass_neighbors
    dstgrid%xlong    = grid%gauss_mass_lon
    dstgrid%xlat     = grid%gauss_mass_lat
    dstgrid%weights  = grid%atmos_mass_2_gauss_mass_weights
    dstgrid%grdnbors = grid%atmos_mass_2_gauss_mass_neighbors

    ! Check local variable and proceed accordingly

    if(is_pinterp) then

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Define local variables

          dstgrid%ncoords      = nlon*nlat
          dstgrid_pres%ncoords = dstgrid%ncoords
          dstgrid_pres%nsig    = nlev

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)
       call mpi_bcast(dstgrid_pres%ncoords,1,mpi_integer,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)
       call mpi_bcast(dstgrid_pres%nsig,1,mpi_integer,mpi_masternode,     &
            & mpi_comm_world,mpi_ierror)
       call mpi_bcast(nlev,1,mpi_integer,mpi_masternode,mpi_comm_world,   &
            & mpi_ierror)

       ! Deallocate memory for local variables

       if(allocated(srcgrid_var)) deallocate(srcgrid_var)
       if(allocated(dstgrid_var)) deallocate(dstgrid_var)

       ! Allocate memory for local variables

       if(.not. allocated(dstgrid%anlysvar))                              &
            & allocate(dstgrid%anlysvar(dstgrid%ncoords,                  &
            & dstgrid%nvertlevs))
       if(.not. allocated(srcgrid_var))                                   &
            & allocate(srcgrid_var(srcgrid%ncoords))
       if(.not. allocated(dstgrid_var))                                   &
            & allocate(dstgrid_var(dstgrid%ncoords))

       ! Initialize local variables
    
       call interpolation_initialize_interppres(srcgrid_pres)
       call interpolation_initialize_interppres(dstgrid_pres)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Compute local variable

          call interpolation_mpas_pressure_profile(srcgrid,dstgrid)

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid%surface_pressure,srcgrid%ncoords,mpi_real,  &
            & mpi_masternode,mpi_comm_world,mpi_ierror)
       call mpi_bcast(srcgrid%pressure_profile,                           &
            & (srcgrid%ncoords*srcgrid%nvertlevs),mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variable

       srcgrid_var = srcgrid%surface_pressure

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)
       
       ! Define local variable

       dstgrid%surface_pressure = dstgrid_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid%surface_pressure,dstgrid%ncoords,mpi_real,  &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Loop through vertical coordinate

       do k = 1, srcgrid%nvertlevs

          ! If on master (root) node (task), define problem and
          ! broadcast variables to each slave (compute) node (task)

          if(mpi_procid .eq. mpi_masternode) then

             ! Define local variable

             srcgrid_var = srcgrid%pressure_profile(k,:)

          end if ! if(mpi_procid .eq. mpi_masternode)

          ! Enable the root task to catch up from I/O and calculations

          call mpi_barrier(mpi_comm_world,mpi_ierror)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,           &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variables

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)
       
          ! If on master (root) node (task), define problem and
          ! broadcast variables to each slave (compute) node (task)

          if(mpi_procid .eq. mpi_masternode) then

             ! Define local variable

             dstgrid%pressure_profile(k,1:dstgrid%ncoords) =               &
                  & dstgrid_var
          
          end if ! if(mpi_procid .eq. mpi_masternode) 

          ! Enable the root task to catch up from I/O and calculations

          call mpi_barrier(mpi_comm_world,mpi_ierror)

       end do ! do k = 1, srcgrid%nvertlevs

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid%pressure_profile,(dstgrid%ncoords*          &
            & srcgrid%nvertlevs),mpi_real,mpi_masternode,mpi_comm_world,  &
            & mpi_ierror)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Loop through vertical coordinate

          do k = 1, nlev

             ! Define local variable

             dstgrid_pres%profile_pres(:,k) = plevs(k)

          end do ! do k = 1, nlev

       end if ! if(mpi_procid .eq. mpi_masternode) 

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid_pres%profile_pres,(dstgrid%ncoords*nlev),   &
            & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Deallocate memory for local variables

       if(allocated(srcgrid_var)) deallocate(srcgrid_var)
       if(allocated(dstgrid_var)) deallocate(dstgrid_var)
       
    end if ! if(is_pinterp)

    ! Check local variable and proceed accordingly

    if(is_zinterp) then

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Define local variables

          dstgrid%ncoords      = nlon*nlat
          dstgrid_hght%ncoords = dstgrid%ncoords
          dstgrid_hght%nsig    = nlev

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)
       call mpi_bcast(dstgrid_hght%ncoords,1,mpi_integer,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)
       call mpi_bcast(dstgrid_hght%nsig,1,mpi_integer,mpi_masternode,     &
            & mpi_comm_world,mpi_ierror)
       call mpi_bcast(nlev,1,mpi_integer,mpi_masternode,mpi_comm_world,   &
            & mpi_ierror)

       ! Deallocate memory for local variables

       if(allocated(srcgrid_var)) deallocate(srcgrid_var)
       if(allocated(dstgrid_var)) deallocate(dstgrid_var)

       ! Allocate memory for local variables

       if(.not. allocated(dstgrid%anlysvar))                              &
            & allocate(dstgrid%anlysvar(dstgrid%ncoords,                  &
            & dstgrid%nvertlevs))
       if(.not. allocated(srcgrid_var))                                   &
            & allocate(srcgrid_var(srcgrid%ncoords))
       if(.not. allocated(dstgrid_var))                                   &
            & allocate(dstgrid_var(dstgrid%ncoords))

       ! Initialize local variables
    
       call interpolation_initialize_interphght(srcgrid_hght)
       call interpolation_initialize_interphght(dstgrid_hght)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Compute local variable

          call interpolation_mpas_height_profile(srcgrid)

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid%surface_height,srcgrid%ncoords,mpi_real,    &
            & mpi_masternode,mpi_comm_world,mpi_ierror)
       call mpi_bcast(srcgrid%height_profile,                             &
            & (srcgrid%ncoords*srcgrid%nvertlevs),mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variable

       srcgrid_var = srcgrid%surface_height

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)
       
       ! Define local variable

       dstgrid%surface_height = dstgrid_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid%surface_height,dstgrid%ncoords,mpi_real,    &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Loop through vertical coordinate

       do k = 1, srcgrid%nvertlevs

          ! If on master (root) node (task), define problem and
          ! broadcast variables to each slave (compute) node (task)

          if(mpi_procid .eq. mpi_masternode) then

             ! Define local variable

             srcgrid_var = srcgrid%height_profile(k,:)

          end if ! if(mpi_procid .eq. mpi_masternode)

          ! Enable the root task to catch up from I/O and calculations

          call mpi_barrier(mpi_comm_world,mpi_ierror)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variables

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)
       
          ! If on master (root) node (task), define problem and
          ! broadcast variables to each slave (compute) node (task)

          if(mpi_procid .eq. mpi_masternode) then

             ! Define local variable

             dstgrid%height_profile(k,1:dstgrid%ncoords) =                 &
                  & dstgrid_var
          
          end if ! if(mpi_procid .eq. mpi_masternode) 

          ! Enable the root task to catch up from I/O and calculations

          call mpi_barrier(mpi_comm_world,mpi_ierror)

       end do ! do k = 1, srcgrid%nvertlevs

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid%height_profile,(dstgrid%ncoords*            &
            & srcgrid%nvertlevs),mpi_real,mpi_masternode,mpi_comm_world,  &
            & mpi_ierror)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Loop through vertical coordinate

          do k = 1, nlev

             ! Define local variable

             dstgrid_hght%profile_hght(:,k) = zlevs(k)

          end do ! do k = 1, nlev

       end if ! if(mpi_procid .eq. mpi_masternode) 

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid_hght%profile_hght,(dstgrid%ncoords*nlev),   &
            & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Deallocate memory for local variables

       if(allocated(srcgrid_var)) deallocate(srcgrid_var)
       if(allocated(dstgrid_var)) deallocate(dstgrid_var)
       
    end if ! if(is_zinterp)

    ! Check local variable and proceed accordingly

    if(is_teinterp) then

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Define local variables

          dstgrid%ncoords      = nlon*nlat
          dstgrid_thee%ncoords = dstgrid%ncoords
          dstgrid_thee%nsig    = nlev

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)
       call mpi_bcast(dstgrid_thee%ncoords,1,mpi_integer,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)
       call mpi_bcast(dstgrid_thee%nsig,1,mpi_integer,mpi_masternode,     &
            & mpi_comm_world,mpi_ierror)
       call mpi_bcast(nlev,1,mpi_integer,mpi_masternode,mpi_comm_world,   &
            & mpi_ierror)

       ! Deallocate memory for local variables

       if(allocated(srcgrid_var)) deallocate(srcgrid_var)
       if(allocated(dstgrid_var)) deallocate(dstgrid_var)

       ! Allocate memory for local variables

       if(.not. allocated(dstgrid%anlysvar))                              &
            & allocate(dstgrid%anlysvar(dstgrid%ncoords,                  &
            & dstgrid%nvertlevs))
       if(.not. allocated(srcgrid_var))                                   &
            & allocate(srcgrid_var(srcgrid%ncoords))
       if(.not. allocated(dstgrid_var))                                   &
            & allocate(dstgrid_var(dstgrid%ncoords))

       ! Initialize local variables
    
       call interpolation_initialize_interpthetae(srcgrid_thee)
       call interpolation_initialize_interpthetae(dstgrid_thee)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Compute local variable

          call interpolation_mpas_thetae_profile(srcgrid)

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid%thetae_profile,                             &
            & (srcgrid%ncoords*srcgrid%nvertlevs),mpi_real,               &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Loop through vertical coordinate

       do k = 1, srcgrid%nvertlevs

          ! If on master (root) node (task), define problem and
          ! broadcast variables to each slave (compute) node (task)

          if(mpi_procid .eq. mpi_masternode) then

             ! Define local variable

             srcgrid_var = srcgrid%thetae_profile(k,:)

          end if ! if(mpi_procid .eq. mpi_masternode)

          ! Enable the root task to catch up from I/O and calculations

          call mpi_barrier(mpi_comm_world,mpi_ierror)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variables

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)
       
          ! If on master (root) node (task), define problem and
          ! broadcast variables to each slave (compute) node (task)

          if(mpi_procid .eq. mpi_masternode) then

             ! Define local variable

             dstgrid%thetae_profile(k,1:dstgrid%ncoords) =                 &
                  & dstgrid_var
          
          end if ! if(mpi_procid .eq. mpi_masternode) 

          ! Enable the root task to catch up from I/O and calculations

          call mpi_barrier(mpi_comm_world,mpi_ierror)

       end do ! do k = 1, srcgrid%nvertlevs

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid%thetae_profile,(dstgrid%ncoords*            &
            & srcgrid%nvertlevs),mpi_real,mpi_masternode,mpi_comm_world,  &
            & mpi_ierror)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Loop through vertical coordinate

          do k = 1, nlev

             ! Define local variable

             dstgrid_thee%profile_thetae(:,k) = televs(k)

          end do ! do k = 1, nlev

       end if ! if(mpi_procid .eq. mpi_masternode) 

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid_thee%profile_thetae,(dstgrid%ncoords*nlev), &
            & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Deallocate memory for local variables

       if(allocated(srcgrid_var)) deallocate(srcgrid_var)
       if(allocated(dstgrid_var)) deallocate(dstgrid_var)
       
    end if ! if(is_teinterp)

    ! Loop through total number of variables

    do l = 1, number_of_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(variable_info_grid(l)%is_profile,1,mpi_logical,     &
            & mpi_masternode,mpi_comm_world,mpi_ierror)
       call mpi_bcast(variable_info_grid(l)%is_derived,1,mpi_logical,     &
            & mpi_masternode,mpi_comm_world,mpi_ierror)
       call mpi_bcast(variable_info_grid(l)%mpas_name,30,mpi_character,   &
            & mpi_masternode,mpi_comm_world,mpi_ierror)
       call mpi_bcast(variable_info_grid(l)%grads_id,10,mpi_character,    &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Allocate memory for local variable

       if(.not. allocated(srcgrid_var))                                   &
            & allocate(srcgrid_var(srcgrid%ncoords))
       if(.not. allocated(dstgrid_var))                                   &
            & allocate(dstgrid_var(dstgrid%ncoords))
       
       ! Check local variable and proceed accordingly

       if(variable_info_grid(l)%is_profile) then

          ! Define local variable
                   
          dstgrid%nvertlevs = srcgrid%nvertlevs

       else   ! if(variable_info_grid(l)%is_profile)

          ! Define local variable
             
          dstgrid%nvertlevs = 1

       end if ! if(variable_info_grid(l)%is_profile)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid%nvertlevs,1,mpi_integer,                    &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Allocate memory for local variable

       if(.not. allocated(workgrid_1d))                                   &
            & allocate(workgrid_1d(srcgrid%ncoords*srcgrid%nvertlevs))
       if(.not. allocated(workgrid_2d))                                   &
            & allocate(workgrid_2d(srcgrid%nvertlevs,srcgrid%ncoords))
       if(.not. allocated(workgrid_2r))                                   &
            & allocate(workgrid_2r(dstgrid%nvertlevs,dstgrid%ncoords))
       if(.not. allocated(workgrid_1r))                                   &
            & allocate(workgrid_1r(dstgrid%ncoords*dstgrid%nvertlevs))

       ! Define local variable accordingly
       
       if(.not. variable_info_grid(l)%is_derived) then

          ! Define local variable
          
          call variable_interface_native(                                 &
               & variable_info_grid(l)%mpas_name,                         &
               & variable_info_grid(l)%is_profile,workgrid_1d)

       else   ! if(.not. variable_info_grid(l)%is_derived)

          ! Define local variable

          call variable_interface_derived(srcgrid,dstgrid,grid,           &
               & variable_info_grid(l)%grads_id,                          &
               & variable_info_grid(l)%is_profile,workgrid_1d,            &
               & workgrid_2r)
             
       end if ! if(.not. variable_info_grid(l)%is_derived)

       ! Check local variable and proceed accordingly
       
       if(variable_info_grid(l)%grads_id .ne. 'psi' .and.                 &
            & variable_info_grid(l)%grads_id .ne. 'chi') then

          ! If on master (root) node (task), define problem and
          ! broadcast variables to each slave (compute) node (task)

          if(mpi_procid .eq. mpi_masternode) then
                
             ! Initialize counting variable

             count = 1

             ! Check local variable and proceed accordingly
          
             if(variable_info_grid(l)%is_profile) then
                
                ! Loop through vertical coordinate
                
                do k = 1, dstgrid%nvertlevs
                   
                   ! Loop through horizontal coordinate
                   
                   do i = 1, srcgrid%ncoords
                      
                      ! Define local variable
                      
                      workgrid_2d(k,i) = workgrid_1d(count)
                      
                      ! Update counting variable
                      
                      count = count + 1
                      
                   end do ! do i = 1, srcgrid%ncoords
                
                end do ! do k = 1, dstgrid%nvertlevs

             else   ! if(variable_info_grid(l)%is_profile)

                ! Loop through horizontal coordinate

                do i = 1, srcgrid%ncoords                    
 
                   ! Define local variable
                   
                   workgrid_2d(1,i) = workgrid_1d(count)
                   
                   ! Update counting variable
                   
                   count = count + 1

                end do ! do i = 1, srcgrid%ncoords

             end if ! if(variable_info_grid(l)%is_profile)

          end if ! if(mpi_procid .eq. mpi_masternode)

          ! Enable the root task to catch up from I/O and calculations

          call mpi_barrier(mpi_comm_world,mpi_ierror)
             
          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(workgrid_2d,(srcgrid%ncoords*dstgrid%nvertlevs),  &
               & mpi_double,mpi_masternode,mpi_comm_world,mpi_ierror)
          
          ! Loop through vertical coordinate
          
          do k = 1, dstgrid%nvertlevs
             
             ! If on master (root) node (task), define problem and
             ! broadcast variables to each slave (compute) node (task)
             
             if(mpi_procid .eq. mpi_masternode) then 
                
                ! Define local variable
                
                srcgrid_var = real(workgrid_2d(k,:))
                
             end if ! if(mpi_procid .eq. mpi_masternode)
             
             ! Enable the root task to catch up from I/O and
             ! calculations
             
             call mpi_barrier(mpi_comm_world,mpi_ierror)
             
             ! Broadcast all necessary variables to slave (compute)
             ! nodes (tasks)
       
             call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,          &
                  & mpi_masternode,mpi_comm_world,mpi_ierror)

             ! Compute local variables

             call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)

             ! If on master (root) node (task), define problem and
             ! broadcast variables to each slave (compute) node (task)
          
             if(mpi_procid .eq. mpi_masternode) then 

                ! Check local variable and proceed accordingly

                if(is_smooth) then

                   ! Compute local variable

                   call interpolation_smoothanalysis(dstgrid_var)

                end if ! if(is_smooth)

                ! Define local variable

                workgrid_2r(k,:) = dstgrid_var(1:dstgrid%ncoords)
             
                ! Check local variable and proceed accordingly

                if(is_pinterp) dstgrid%anlysvar(1:dstgrid%ncoords,k) =     &
                     & dstgrid_var(1:dstgrid%ncoords)   
                if(is_zinterp) dstgrid%anlysvar(1:dstgrid%ncoords,k) =     &
                     & dstgrid_var(1:dstgrid%ncoords)             
                if(is_teinterp) dstgrid%anlysvar(1:dstgrid%ncoords,k) =    &
                     & dstgrid_var(1:dstgrid%ncoords)         

             end if ! if(mpi_procid .eq. mpi_masternode)

             ! Enable the root task to catch up from I/O and
             ! calculations
             
             call mpi_barrier(mpi_comm_world,mpi_ierror)

          end do ! do k = 1, dstgrid%nvertlevs

       else   ! if(variable_info_grid(l)%grads_id .ne. 'psi' .and.         &
              ! variable_info_grid(l)%grads_id .ne. 'chi')

          ! Loop through vertical coordinate
          
          do k = 1, dstgrid%nvertlevs

             ! If on master (root) node (task), define problem and
             ! broadcast variables to each slave (compute) node (task)
             
             if(mpi_procid .eq. mpi_masternode) then 

                ! Check local variable and proceed accordingly
          
                if(is_pinterp) dstgrid%anlysvar(1:dstgrid%ncoords,k) =     &
                     & workgrid_2r(k,1:dstgrid%ncoords)
                if(is_zinterp) dstgrid%anlysvar(1:dstgrid%ncoords,k) =     &
                     & workgrid_2r(k,1:dstgrid%ncoords)
                if(is_teinterp) dstgrid%anlysvar(1:dstgrid%ncoords,k) =    &
                     & workgrid_2r(k,1:dstgrid%ncoords)

             end if ! if(mpi_procid .eq. mpi_masternode)

             ! Enable the root task to catch up from I/O and
             ! calculations
             
             call mpi_barrier(mpi_comm_world,mpi_ierror)

          end do ! do k = 1, dstgrid%nvertlevs

       end if ! if(variable_info_grid(l)%grads_id .ne. 'psi' .and.         &
              ! variable_info_grid(l)%grads_id .ne. 'chi')

       ! Check local variable and proceed accordingly

       if(is_pinterp .or. is_zinterp .or. is_teinterp) then

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(dstgrid%anlysvar,                                 &
               & (dstgrid%ncoords*dstgrid%nvertlevs),mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end if ! if(is_pinterp .or. is_zinterp .or. is_teinterp)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Allocate memory for local variable accordingly

          if((is_pinterp .or. is_zinterp .or. is_teinterp) .and.           &
               & variable_info_grid(l)%is_profile) then
             
             ! Allocate memory for local variable
             
             if(.not. allocated(variable_info_grid(l)%array))              &
                  & allocate(variable_info_grid(l)%array(                  &
                  & dstgrid%ncoords*nlev))

          else if(.not. is_pinterp .and. .not. is_zinterp .and. .not.      &
               & is_teinterp .and. variable_info_grid(l)%is_profile) then

             ! Allocate memory for local variable
                   
             if(.not. allocated(variable_info_grid(l)%array))              &
                  & allocate(variable_info_grid(l)%array(                  &
                  & dstgrid%ncoords*dstgrid%nvertlevs))

          else   ! if(is_pinterp .or. is_zinterp .or. is_teinterp)

             ! Allocate memory for local variable
                   
             if(.not. allocated(variable_info_grid(l)%array))              &
                  & allocate(variable_info_grid(l)%array(                  &
                  & dstgrid%ncoords*dstgrid%nvertlevs))

          end if ! if(is_pinterp .or. is_zinterp .or. is_teinterp)

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Check local variable and proceed accordingly

       if(is_pinterp .and. variable_info_grid(l)%is_profile) then

          ! Compute local variable

          call interpolation_interppres_mpi(dstgrid_pres,dstgrid)
       
       end if ! if(is_pinterp .and. variable_info_grid(l)%is_profile)

       ! Check local variable and proceed accordingly

       if(is_zinterp .and. variable_info_grid(l)%is_profile) then

          ! Compute local variable

          call interpolation_interphght_mpi(dstgrid_hght,dstgrid)
       
       end if ! if(is_zinterp .and. variable_info_grid(l)%is_profile)

       ! Check local variable and proceed accordingly

       if(is_teinterp .and. variable_info_grid(l)%is_profile) then

          ! Compute local variable

          call interpolation_interpthetae_mpi(dstgrid_thee,dstgrid)
       
       end if ! if(is_teinterp .and. variable_info_grid(l)%is_profile)

       ! If on master (root) node (task), define problem and
       ! broadcast variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Initialize local variable

          array_start = 1
          array_end   = array_start + (dstgrid%ncoords - 1)

          ! Check local variable and proceed accordingly

          if(.not. is_pinterp .and. .not. is_zinterp .and. .not.          &
               & is_teinterp) then

             ! Loop through vertical coordinate

             do k = 1, dstgrid%nvertlevs                    

                ! Define local variable

                variable_info_grid(l)%array(array_start:array_end) =      &
                     & workgrid_2r(k,:)

                ! Update local variables

                array_start = array_end + 1
                array_end   = array_start + (dstgrid%ncoords - 1)

             end do ! do k = 1, dstgrid%nvertlevs

          else if(.not. variable_info_grid(l)%is_profile) then ! if(.not. is_pinterp
                                                               ! .and. .not. is_zinterp
                                                               ! .and. .not. is_teinterp)

             ! Define local variable

             variable_info_grid(l)%array(array_start:array_end) =            &       
                  & workgrid_2r(1,:)             

          else if(is_pinterp) then                             ! if(.not. is_pinterp
                                                               ! .and. .not. is_zinterp
                                                               ! .and. .not. is_teinterp)
             
             ! Loop through vertical coordinate
             
             do k = 1, nlev
                
                ! Define local variable
                
                variable_info_grid(l)%array(array_start:array_end) =         &
                     & dstgrid_pres%profile_var(:,k)
                
                ! Update local variables

                array_start = array_end + 1
                array_end   = array_start + (dstgrid%ncoords - 1)
                
             end do ! do k = 1, nlev
             
          else if(is_zinterp) then                             ! if(.not. is_pinterp
                                                               ! .and. .not. is_zinterp
                                                               ! .and. .not. is_teinterp)

             ! Loop through vertical coordinate
             
             do k = 1, nlev
                
                ! Define local variable
                
                variable_info_grid(l)%array(array_start:array_end) =         &
                     & dstgrid_hght%profile_var(:,k)
                
                ! Update local variables

                array_start = array_end + 1
                array_end   = array_start + (dstgrid%ncoords - 1)
                
             end do ! do k = 1, nlev

          else if(is_teinterp) then                            ! if(.not. is_pinterp
                                                               ! .and. .not. is_zinterp
                                                               ! .and. .not. is_teinterp)

             ! Loop through vertical coordinate
             
             do k = 1, nlev
                
                ! Define local variable
                
                variable_info_grid(l)%array(array_start:array_end) =         &
                     & dstgrid_thee%profile_var(:,k)
                
                ! Update local variables

                array_start = array_end + 1
                array_end   = array_start + (dstgrid%ncoords - 1)
                
             end do ! do k = 1, nlev
             
          end if                                               ! if(.not. is_pinterp
                                                               ! .and. .not. is_zinterp
                                                               ! .and. .not. is_teinterp)

          ! Check local variable and proceed accordingly

          if(variable_info_grid(l)%grads_id .eq. 'qv' .or.                 &
               & variable_info_grid(l)%grads_id .eq. 'rh') then

             ! Compute local variable

             call calculate_clipvariable((nlon*nlat),                      &
                  & variable_info_grid(l)%zdim,                            &
                  & variable_info_grid(l)%array)

          end if ! if(variable_info_grid(l)%grads_id .eq. 'qv' .or.        &
                 ! variable_info_grid(l)%grads_id .eq. 'rh')

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Deallocate memory for local variables

       if(allocated(srcgrid_var)) deallocate(srcgrid_var)
       if(allocated(dstgrid_var)) deallocate(dstgrid_var)
       if(allocated(workgrid_1d)) deallocate(workgrid_1d)
       if(allocated(workgrid_2d)) deallocate(workgrid_2d)
       if(allocated(workgrid_1r)) deallocate(workgrid_1r)
       if(allocated(workgrid_2r)) deallocate(workgrid_2r)
       
    end do ! do l = 1, number_of_var

    ! Deallocate memory for local variables

    call interpolation_cleanup_interppres(srcgrid_pres)
    call interpolation_cleanup_interppres(dstgrid_pres)
    call interpolation_cleanup_interphght(srcgrid_hght)
    call interpolation_cleanup_interphght(dstgrid_hght)
    call interpolation_cleanup_interpthetae(srcgrid_thee)
    call interpolation_cleanup_interpthetae(dstgrid_thee)
    call interpolation_cleanup_grid(srcgrid)
    call interpolation_cleanup_grid(dstgrid)
    call interpolation_cleanup_task_balance(srcgrid)
    call interpolation_cleanup_task_balance(dstgrid)
    if(allocated(dstgrid%anlysvar)) deallocate(dstgrid%anlysvar)

    !=====================================================================

  end subroutine variable_interface_process

  !=======================================================================

  ! variable_interface_native.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_native(mpas_name,is_profile,src_variable)

    ! Define variables passed to routine

    character(len=30)                                             :: mpas_name
    logical                                                       :: is_profile

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_nlevels*mpas_ncoords) :: src_variable

    ! Define variables computed within routine

    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: workgrid_2d
    real(r_double),          dimension(mpas_ncoords)              :: workgrid_1d

    ! Define counting variables

    integer                                                       :: i, j, k, l
    integer                                                       :: count

    !=====================================================================

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),                     &
         & mode=nf90_nowrite,ncid=ncfileid)

    ! Check local variable and proceed accordingly

    if(.not. is_profile) then
          
       ! Define local variables

       ncstatus = nf90_inq_varid(ncfileid,                                 &
            & trim(adjustl(mpas_name)),ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid_1d)
       
    else   ! if(.not. is_profile)

       ! Define local variables
       
       ncstatus = nf90_inq_varid(ncfileid,                                 &
            & trim(adjustl(mpas_name)),ncvarid)
       ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid_2d)

    end if ! if(.not. is_profile)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize counting variable
       
    count = 1
       
    ! Check local variable and proceed accordingly
       
    if(.not. is_profile) then
          
       ! Loop through horizontal coordinate
          
       do i = 1, mpas_ncoords
             
          ! Define local variable
             
          src_variable(count) = workgrid_1d(i)
             
          ! Update counting variable
             
          count = count + 1
             
       end do ! do i = 1, mpas_ncoords
          
    else   ! if(.not. is_profile)
          
       ! Loop through vertical coordinate
          
       do k = 1, mpas_nlevels
             
          ! Loop through horizontal coordinate
             
          do i = 1, mpas_ncoords
                
             ! Define local variable
                
             src_variable(count) = workgrid_2d(k,i)
                
             ! Update counting variable
             
             count = count + 1
             
          end do ! do i = 1, mpas_ncoords
          
       end do ! do k = 1, mpas_nlevels
       
    end if ! if(.not. is_profile)

    !=====================================================================

  end subroutine variable_interface_native

  !=======================================================================

  ! variable_interface_derived.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_derived(srcgrid,dstgrid,grid,grads_id,     &
       & is_profile,src_variable,dst_variable)

    ! Define variables passed to routine

    type(interpgrid)                                                      :: srcgrid
    type(interpgrid)                                                      :: dstgrid
    type(grid_interface)                                                  :: grid
    character(len=10)                                                     :: grads_id
    logical                                                               :: is_profile

    ! Define variables returned by routine

    real(r_double),          dimension(srcgrid%nvertlevs*srcgrid%ncoords) :: src_variable
    real(r_kind),            dimension(dstgrid%nvertlevs,dstgrid%ncoords) :: dst_variable

    !=====================================================================

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grads_id,10,mpi_character,mpi_masternode,               &
         & mpi_comm_world,mpi_ierror)

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Compute local variable accordingly

       if(grads_id .eq. 'pslp') then
          
          ! Compute local variable
          
          call calculate_sealevelpressure(src_variable)
          
       end if ! if(grads_id .eq. 'pslp')

       ! Compute local variable accordingly

       if(grads_id .eq. 'terrain') then
          
          ! Compute local variable
          
          call calculate_terrain(src_variable)
          
       end if ! if(grads_id .eq. 'ter')

       ! Compute local variable accordingly

       if(grads_id .eq. 'cori') then
          
          ! Compute local variable
          
          call calculate_coriolis(src_variable)
          
       end if ! if(grads_id .eq. 'cori')
       
       ! Compute local variable accordingly
       
       if(grads_id .eq. 'vort') then
          
          ! Compute local variable
          
          call calculate_vorticity(src_variable)
          
       end if ! if(grads_id .eq. 'vort')

       ! Compute local variable accordingly
       
       if(grads_id .eq. 'divg') then
          
          ! Compute local variable
          
          call calculate_divergence(src_variable)
          
       end if ! if(grads_id .eq. 'divg')
       
       ! Compute local variable accordingly
       
       if(grads_id .eq. 'rh') then
          
          ! Compute local variable
          
          call calculate_relativehumidity(src_variable)
          
       end if ! if(grads_id .eq. 'rh')
       
       ! Compute local variable accordingly
       
       if(grads_id .eq. 't') then
          
          ! Compute local variable
          
          call calculate_temperature(src_variable)
          
       end if ! if(grads_id .eq. 't')

       ! Compute local variable accordingly
       
       if(grads_id .eq. 'tv') then
          
          ! Compute local variable
          
          call calculate_virtualtemperature(src_variable)
          
       end if ! if(grads_id .eq. 'tv')
       
       ! Compute local variable accordingly

       if(grads_id .eq. 'thetae') then
          
          ! Compute local variable
          
          call calculate_thetae(src_variable)
          
       end if ! if(grads_id .eq. 'thetae')

      ! Compute local variable accordingly

       if(grads_id .eq. 'pv') then
          
          ! Compute local variable
          
          call calculate_potentialvorticity(src_variable)
          
       end if ! if(grads_id .eq. 'pv')

      ! Compute local variable accordingly

       if(grads_id .eq. 'dbz') then
          
          ! Compute local variable
          
          call calculate_reflectivity_dbz(src_variable)
          
       end if ! if(grads_id .eq. 'dbz')

      ! Compute local variable accordingly

       if(grads_id .eq. 'cdbz') then
          
          ! Compute local variable
          
          call calculate_reflectivity_composite(src_variable)
          
       end if ! if(grads_id .eq. 'cdbz')

      ! Compute local variable accordingly

       if(grads_id .eq. 'ght') then
          
          ! Compute local variable
          
          call calculate_geopotential(src_variable)
          
       end if ! if(grads_id .eq. 'ght')

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Compute local variable accordingly

    if(grads_id .eq. 'psi') then

       ! Compute local variable

       call calculate_streamfunction(srcgrid,dstgrid,dst_variable)
    
    end if ! if(grads_id .eq. 'psi') 

    ! Compute local variable accordingly

    if(grads_id .eq. 'chi') then

       ! Compute local variable

       call calculate_velocitypotential(srcgrid,dstgrid,dst_variable)
    
    end if ! if(grads_id .eq. 'chi') 

    !=====================================================================

  end subroutine variable_interface_derived

  !=======================================================================

  ! calculate_geopotential.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_geopotential(src_variable)

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_nlevels*mpas_ncoords)   :: src_variable

    ! Define variables computed within routine

    real(r_double),          dimension(mpas_nlevels+1,mpas_ncoords) :: height
    real(r_double),          dimension(mpas_ncoords)                :: xlat
    real(r_double),          dimension(mpas_ncoords)                :: grav_lat

    ! Define counting variables

    integer                                                         :: i, j, k
    integer                                                         :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)
    src_variable = 1.e30

    ! Print message to user

    if(debug) write(6,500) 

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'zgrid',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,height) 
    ncstatus = nf90_inq_varid(ncfileid,'latCell',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,xlat) 

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize counting variable

    count = 1

    ! Loop through local variable

    do k = 1, mpas_nlevels

       ! Loop through local variable

       do i = 1, mpas_ncoords

          ! Compute local variables

          grav_lat(i)         = 9.780327*(1.0 + 0.053024*(sin(xlat(i))*     &
               & sin(xlat(i))) - 0.0000025*(sin(2.0*xlat(i))*               &
               & sin(2.0*xlat(i))))
          src_variable(count) = grav_lat(i)*((height(k,i) +                 &
               & height(k+1,i))/2.0)
          
          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

    ! Define format statements

500 format('CALCULATE_GEOPOTENTIAL: Computing MPAS geopotential.')

    !=====================================================================

  end subroutine calculate_geopotential

  !=======================================================================

  ! calculate_terrain.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_terrain(src_variable)

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_ncoords)                :: src_variable

    ! Define variables computed within routine

    real(r_double),          dimension(mpas_nlevels+1,mpas_ncoords) :: height

    ! Define counting variables

    integer                                                         :: i, j, k

    !=====================================================================

    ! Initialize local variable

    src_variable = 1.e30

    ! Print message to user

    if(debug) write(6,500) 

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'zgrid',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,height) 

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Define local variable

    src_variable(1:mpas_ncoords) = height(1,1:mpas_ncoords)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

    ! Define format statements

500 format('CALCULATE_TERRAIN: Computing MPAS terrain elevation.')

    !=====================================================================

  end subroutine calculate_terrain

  !=======================================================================

  ! calculate_sealevelpressure.f90: This subroutine will compute the 
  ! pressure reduced to sea-level; the subroutine proceeds as follows:
  !
  ! (1) Find least zeta level that is PCONST (e.g., 10000) Pa above
  !     the surface to extrapolate a surface pressure and temperature,
  !     which is supposed to reduce the effect of the diurnal heating
  !     cycle in the pressure field
  !
  ! (2) Get temperature PCONST Pa above surface to extrapolate the
  !     temperature at the surface and down to sea level
  !
  ! (3) Compute a correction to the sea level temperature if both the
  !     surface and sea level temperatures are *too* hot
  !
  ! (4) Compute the pressure reduced to sea-level
  !
  !-----------------------------------------------------------------------

  subroutine calculate_sealevelpressure(src_variable)

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_ncoords)                :: src_variable

    ! Define variables computed within routine

    logical                                                         :: found
    logical                                                         :: l1
    logical                                                         :: l2
    logical                                                         :: l3
    real(r_double),          dimension(mpas_nlevels+1,mpas_ncoords) :: height
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords)   :: pressure
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords)   :: temperature
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords)   :: qvapor
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords)   :: workgrid
    real(r_double),          dimension(mpas_ncoords)                :: ter
    real(r_double),          dimension(mpas_ncoords)                :: psfc
    real(r_double),          dimension(mpas_ncoords)                :: t_surf
    real(r_double),          dimension(mpas_ncoords)                :: t_sealevel
    real(r_double),          dimension(mpas_nlevels)                :: rdzw
    real(r_double)                                                  :: p_at_pconst
    real(r_double)                                                  :: t_at_pconst
    real(r_double)                                                  :: z_at_pconst
    real(r_double)                                                  :: plo
    real(r_double)                                                  :: phi
    real(r_double)                                                  :: tlo
    real(r_double)                                                  :: thi
    real(r_double)                                                  :: zlo
    real(r_double)                                                  :: zhi
    integer,                 dimension(mpas_ncoords)                :: level
    integer                                                         :: klo
    integer                                                         :: khi

    ! Define counting variables

    integer                                                         :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)
    src_variable = 1.e30

    ! Print message to user

    if(debug) write(6,500) 

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = pressure + workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'theta',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Compute local variable

    temperature = workgrid/((100000.0/pressure)**rd_over_cp_mass)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'zgrid',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,height) 

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'qv',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,qvapor)

   ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'rdzw',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,rdzw)

   ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'surface_pressure',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,psfc)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Loop through all horizontal coordinates and proceed accordingly

    do i = 1, mpas_ncoords

       ! Define local variable

       level(i) = -1
       k        = 1
       found    = .false.

       ! Define local variable

       ter(i) = height(1,i)

       ! Check local variable and proceed accordingly

       do while((.not. found) .and. (k .le. mpas_nlevels))

          ! Check local variable and proceed accordingly

          if(pressure(k,i) .lt. pressure(k,i) - 10000.0) then

             ! Define local variables

             level(i) = k
             found    = .true.

          end if ! if(pressure(k,i) .lt. pressure(k,i) - 10000.0)

          ! Update local variable

          k = k + 1

       end do ! do while((.not. found) .and. (k .le. mpas_nlevels))

    end do ! do i = 1, mpas_ncoords

    ! Loop through all horizontal coordinates and proceed accordingly

    do i = 1, mpas_ncoords

       ! Define local variables

       klo = max(level(i)-1,1                )
       khi = min(klo + 1   , mpas_nlevels - 1)
       plo = pressure(klo,i)
       phi = pressure(khi,i)
       tlo = temperature(klo,i)*(1.0 + 0.608*qvapor(klo,i))
       thi = temperature(khi,i)*(1.0 + 0.608*qvapor(khi,i))
       zlo = height(klo,i)
       zhi = height(khi,i)

       ! Compute local variables

       p_at_pconst   = psfc(i) - 10000.0
       t_at_pconst   = thi - (thi - tlo)*log(p_at_pconst/phi)*log(plo/phi)
       z_at_pconst   = zhi - (zhi - zlo)*log(p_at_pconst/phi)*log(plo/phi)
       t_surf(i)     = t_at_pconst*(psfc(i)/10000.0)**(gamma*rd/grav)
       t_sealevel(i) = t_at_pconst+gamma*z_at_pconst

       ! Define local variables

       l1 = t_sealevel(i) .lt. (273.16+17.5)
       l2 = t_surf(i)     .le. (273.16+17.5)
       l3 = .not. l1

       ! Check local variables and proceed accordingly

       if(l2 .and. l3) then

          ! Define local variable

          t_sealevel(i) = (273.16+17.5)

       else   ! if(l2 .and. l3)

          ! Define local variable

          t_sealevel(i) = (273.16+17.5) - 0.005*(t_surf(i) -               &
               & (273.16+17.5))**2.0

       end if ! if(l2 .and. l3)

       ! Compute local variable

       src_variable(i) = psfc(i)*exp((2.0*grav*ter(i))/(rd*(t_sealevel(i)  &
            & +t_surf(i))))

    end do ! do i = 1, mpas_ncoords

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

    ! Define format statements

500 format('CALCULATE_SEALEVELPRESSURE: Computing MPAS sea-level ',        &
         & 'pressure.')

    !=====================================================================

  end subroutine calculate_sealevelpressure

  !=======================================================================

  ! calculate_landseaicemask.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_landseaicemask(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_ncoords)                :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(mpas_ncoords)                :: seaice
    integer,                   dimension(mpas_ncoords)                :: landmask

    !=====================================================================

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'seaice',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,seaice)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'landmask',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,landmask)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Compute local variable

    src_variable = dble(landmask)

    ! Rescale local variable accordingly

    where(seaice .eq. 1.0) src_variable = 2.0

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_landseaicemask

  !=======================================================================

  ! calculate_vorticity.f90: This subroutine will compute the relative
  ! vorticity profile from the MPAS native vorticity which is defined
  ! along the cell vertices; the 'reconstructed' relative vorticity is
  ! the mean value from the vorticity values defined on the respective
  ! cell's vertices

  !-----------------------------------------------------------------------

  subroutine calculate_vorticity(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)               :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(:,:),                        allocatable :: vorticity
    real(r_double),            dimension(:,:),                        allocatable :: kiteAreasOnVertex
    real(r_double),            dimension(:,:),                        allocatable :: src_vort
    real(r_double),            dimension(:),                          allocatable :: areaCell
    integer,                   dimension(:,:),                        allocatable :: cellsOnVertex
    integer                                                                       :: nVertices
    integer                                                                       :: vertexDegree
    integer                                                                       :: iCell

    ! Define counting variables

    integer                                                                       :: i, j, k, l
    integer                                                                       :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    ! Print message to user

    if(debug) write(6,500) 

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define local variables

    ncstatus = nf90_inq_dimid(ncfileid,'nVertices',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,len=nVertices)
    ncstatus = nf90_inq_dimid(ncfileid,'vertexDegree',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,len=vertexDegree)

    ! Allocate memory for local variables

    if(.not. allocated(vorticity))                                         &
         & allocate(vorticity(mpas_nlevels,nVertices))
    if(.not. allocated(src_vort))                                          &
         & allocate(src_vort(mpas_nlevels,mpas_ncoords))
    if(.not. allocated(kiteAreasOnVertex))                                 &
         & allocate(kiteAreasOnVertex(vertexDegree,nVertices))
    if(.not. allocated(cellsOnVertex))                                     &
         & allocate(cellsOnVertex(vertexDegree,nVertices))
    if(.not. allocated(areaCell))                                          &
         & allocate(areaCell(mpas_ncoords))

    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'vorticity',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,vorticity)
    ncstatus = nf90_inq_varid(ncfileid,'kiteAreasOnVertex',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,kiteAreasOnVertex)
    ncstatus = nf90_inq_varid(ncfileid,'cellsOnVertex',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,cellsOnVertex)
    ncstatus = nf90_inq_varid(ncfileid,'areaCell',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,areaCell)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize local variables

    src_vort     = 0.0
    src_variable = 1.e30

    ! Loop through total number of vertices

    do j = 1, nVertices

       ! Loop thought total number of vertex degrees

       do i = 1, vertexDegree

          ! Define local variable

          iCell = cellsOnVertex(i,j)

          ! Loop through vertical coordinate

          do k = 1, mpas_nlevels

             ! Compute local variable

             src_vort(k,iCell) = src_vort(k,iCell) +                       &
                  & kiteAreasOnVertex(i,j)*vorticity(k,j)/areaCell(iCell)

          end do ! do k = 1, mpas_nlevels

       end do ! do i = 1, vertexDegree

    end do ! do j = 1, nVertices

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, mpas_nlevels

       ! Loop through horizontal coordinate

       do i = 1, mpas_ncoords

          ! Define local variable

          src_variable(count) = src_vort(k,i)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(vorticity))         deallocate(vorticity)
    if(allocated(src_vort))          deallocate(src_vort)
    if(allocated(cellsOnVertex))     deallocate(cellsOnVertex)
    if(allocated(kiteAreasOnVertex)) deallocate(kiteAreasOnVertex)
    if(allocated(areaCell))          deallocate(areaCell)
       
    !=====================================================================

    ! Define format statements

500 format('CALCULATE_VORTICITY: Computing MPAS relative vorticity.')

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_vorticity

  !=======================================================================

  ! calculate_divergence.f90: This subroutine will ingest the divergence
  ! profile on the native MPAS grid

  !-----------------------------------------------------------------------

  subroutine calculate_divergence(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)               :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(:,:),                        allocatable :: divergence

    ! Define counting variables

    integer                                                                       :: i, j, k, l
    integer                                                                       :: count

    !=====================================================================

    ! Print message to user

    if(debug) write(6,500) 

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Allocate memory for local variables

    if(.not. allocated(divergence))                                        &
         & allocate(divergence(mpas_nlevels,mpas_ncoords))

    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'divergence',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,divergence)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize local variable

    src_variable = 1.e30

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, mpas_nlevels

       ! Loop through horizontal coordinate

       do i = 1, mpas_ncoords

          ! Define local variable

          src_variable(count) = divergence(k,i)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(divergence)) deallocate(divergence)
       
    !=====================================================================

    ! Define format statements

500 format('CALCULATE_DIVERGENCE: Computing MPAS divergence.')

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_divergence

  !=======================================================================

  ! calculate_relativehumidity.f90: This subroutine will compute the
  ! relative humidity profile from the model native potential
  ! temperature, water vapor, and pressure profiles

  !-----------------------------------------------------------------------

  subroutine calculate_relativehumidity(src_variable)

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_nlevels*mpas_ncoords) :: src_variable

    ! Define variables computed within routine

    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: pressure
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: temperature
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: qvapor
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: workgrid
    real(r_double)                                                :: es
    real(r_double)                                                :: qs

    ! Define counting variables

    integer                                                       :: i, j, k
    integer                                                       :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)
    src_variable = 1.e30    

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = pressure + workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'qv',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    qvapor = workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'theta',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Compute local variable

    temperature = workgrid/((100000.0/pressure)**rd_over_cp_mass)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1

    ! Loop through all vertical coordinates and proceed accordingly

    do k = 1, mpas_nlevels

       ! Loop through all horizontal coordinates and proceed
       ! accordingly

       do i = 1, mpas_ncoords

          ! Check local variable and compute local variable
          ! accordingly

          if(temperature(k,i) .le. 273.15) then

             ! Compute local variable

             es = dble(6.11)*exp(dble(22.514)-(dble(6150.0)/               &
                  & temperature(k,i)))

          else   ! if(temperature(k,i) .le. 273.15)

             ! Compute local variable

             es = dble(6.112)*exp(dble(17.67)*((temperature(k,i) -         &
                  & dble(273.15))/(temperature(k,i) - dble(29.65))))
             
          end if ! if(temperature(k,i) .le. 273.15)

          ! Compute local variables

          qs                  = dble(0.622)*es/((pressure(k,i)/            &
               & dble(100.0)) - es)
          src_variable(count) = dble(100.0)*(qvapor(k,i)/qs)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    ! Rescale local variable accordingly

    where(src_variable .gt. dble(100.0)) src_variable = dble(100.0)
    where(src_variable .lt. dble(0.0  )) src_variable = dble(0.0  )

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_relativehumidity

  !=======================================================================

  ! calculate_temperature.f90: This subroutine will compute the
  ! temperature profile from the model native potential temperature
  ! and pressure profiles

  !-----------------------------------------------------------------------

  subroutine calculate_temperature(src_variable)

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_nlevels*mpas_ncoords) :: src_variable

    ! Define variables computed within routine

    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: pressure
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: workgrid

    ! Define counting variables

    integer                                                       :: i, j, k
    integer                                                       :: count
    
    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)
    src_variable = 1.e30

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = pressure + workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'theta',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1

    ! Loop through all vertical coordinates and proceed accordingly

    do k = 1, mpas_nlevels

       ! Loop through all horizontal coordinates and proceed
       ! accordingly

       do i = 1, mpas_ncoords

          ! Compute local variable

          src_variable(count) = workgrid(k,i)/((100000.0/pressure(k,i))    &
               & **rd_over_cp_mass)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_temperature

  !=======================================================================

  ! calculate_virtualtemperature.f90: This subroutine will compute the
  ! virtual temperature profile from the model native potential
  ! temperature, moisture, and pressure profiles

  !-----------------------------------------------------------------------

  subroutine calculate_virtualtemperature(src_variable)

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_nlevels*mpas_ncoords) :: src_variable

    ! Define variables computed within routine

    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: pressure
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: theta
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: workgrid

    ! Define counting variables

    integer                                                       :: i, j, k
    integer                                                       :: count
    
    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)
    src_variable = 1.e30

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = pressure + workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'theta',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    theta = workgrid

   ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'qv',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1

    ! Loop through all vertical coordinates and proceed accordingly

    do k = 1, mpas_nlevels

       ! Loop through all horizontal coordinates and proceed
       ! accordingly

       do i = 1, mpas_ncoords

          ! Compute local variables

          src_variable(count) = theta(k,i)/((100000.0/pressure(k,i))      &
               & **rd_over_cp_mass)
          src_variable(count) = src_variable(count)*(1.0 +                &
               & 0.61*workgrid(k,i))

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_virtualtemperature

  !=======================================================================

  ! calculate_thetae.f90: This subroutine will compute the equivalent
  ! potential temperature profile from the model native potential
  ! temperature and pressure profiles

  !-----------------------------------------------------------------------

  subroutine calculate_thetae(src_variable)

    ! Define variable returned by routine

    real(r_double),          dimension(mpas_nlevels*mpas_ncoords) :: src_variable

    ! Define variables computed within routine

    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: pressure
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: theta
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: qvapor
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: temperature
    real(r_double),          dimension(mpas_nlevels,mpas_ncoords) :: workgrid
    real(r_double)                                                :: vappres
    real(r_double)                                                :: tlcl

    ! Define counting variables

    integer                                                       :: i, j, k
    integer                                                       :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)
    src_variable = 1.e30

    !---------------------------------------------------------------------

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)

    ! Define local variable

    pressure = pressure + workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'theta',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Compute local variable

    theta = workgrid

    ! Define and ingest local variables

    ncstatus = nf90_inq_varid(ncfileid,'qv',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Define local variable

    qvapor = workgrid

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    ! Compute local variable

    temperature = theta/((100000.0/pressure)**rd_over_cp_mass)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1

    ! Loop through all vertical coordinates and proceed accordingly

    do k = 1, mpas_nlevels

       ! Loop through all horizontal coordinates and proceed
       ! accordingly

       do i = 1, mpas_ncoords

          ! Compute local variables

          vappres = pressure(k,i)*qvapor(k,i)/(dble(eps)+qvapor(k,i))/     &
               & dble(100.0)
          tlcl    = 55.0 + 2840.0/(3.5*log(temperature(k,i)) -             &
               & log(vappres) - 4.805)
          
          ! Rescale local variable accordingly

          tlcl = min(tlcl,temperature(k,i))

          ! Compute local variable

          src_variable(count) = theta(k,i)*exp(((3376./tlcl)-2.54)*        &
               & qvapor(k,i)*(1.0 + 0.81*qvapor(k,i)))

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_thetae

  !=======================================================================

  ! calculate_streamfunction.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_streamfunction(srcgrid,dstgrid,dst_variable)

    ! Define variables passed to routine

    type(interpgrid)                                                                    :: srcgrid
    type(interpgrid)                                                                    :: dstgrid
    type(grid_interface)                                                                :: grid

    ! Define variable returned by routine

    real(r_kind),              dimension(dstgrid%nvertlevs,dstgrid%ncoords)             :: dst_variable

    ! Define variables computed within routine

    real(r_double),            dimension(srcgrid%nvertlevs*srcgrid%ncoords)             :: workgrid
    real(r_kind),              dimension(nlon,nlat,dstgrid%nvertlevs)                   :: dst_vort
    real(r_kind),              dimension(nlon,nlat,dstgrid%nvertlevs)                   :: dst_psi
    real(r_kind),              dimension(nlon+1,nlat+1)                                 :: rhs
    real(r_kind),              dimension(nlon,nlat)                                     :: vort
    real(r_kind),              dimension(nlon,nlat)                                     :: psi
    real(r_kind),              dimension(nlat+1)                                        :: xbdl 
    real(r_kind),              dimension(nlat+1)                                        :: xbdr 
    real(r_kind),              dimension(nlon+1)                                        :: ybdt 
    real(r_kind),              dimension(nlon+1)                                        :: ybdb 
    real(r_kind),              dimension(10000000)                                      :: work
    real(r_kind)                                                                        :: perturb
    real(r_kind)                                                                        :: xmin
    real(r_kind)                                                                        :: xmax
    real(r_kind)                                                                        :: ymin
    real(r_kind)                                                                        :: ymax
    integer                                                                             :: ierror

    ! Define counting variables

    integer                                                                             :: i, j, k, l
    integer                                                                             :: count
    integer                                                                             :: dst_count

    !=====================================================================

    ! Initialize local variable

    dst_variable = 1.e30

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize local variables
       
       call init_constants_derived()
       call init_constants(.true.)
       perturb = 0.0

       ! Compute local variable

       call calculate_vorticity(workgrid)

       ! Compute local variables

       xmin = 1.0
       xmax = dx*cos((sum(dstgrid%xlat)/real(dstgrid%ncoords)))*           &
            & ((2.0*pi*rearth_equator)/360.0)
       ymin = 1.0
       ymax = dy*cos((sum(dstgrid%xlat)/real(dstgrid%ncoords)))*111321.0

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(workgrid,(mpas_ncoords*mpas_nlevels),mpi_real,          &
         & mpi_masternode,mpi_comm_world,mpi_ierror)       
    call mpi_bcast(xmin,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(xmax,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(ymin,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(ymax,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(perturb,1,mpi_real,mpi_masternode,mpi_comm_world,       &
         & mpi_ierror)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, srcgrid%nvertlevs

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then
       
          ! Loop through horizontal coordinate
                   
          do i = 1, srcgrid%ncoords
                      
             ! Define local variable
             
             srcgrid_var(i) = real(workgrid(count))
             
             ! Update counting variable
             
             count = count + 1
             
          end do ! do i = 1, srcgrid%ncoords

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)

       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Initialize local variable

          dst_count = 1

          ! Loop through meridional horizontal coordinate

          do j = 1, nlat
       
             ! Loop through zonal horizontal coordinate
                   
             do i = 1, nlon
                      
                ! Define local variable
             
                dst_vort(i,j,k) = dstgrid_var(dst_count)
             
                ! Update counting variable
             
                dst_count = dst_count + 1
             
             end do ! do i = 1, nlon

          end do !  do j = 1, nlat

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

    end do ! do k = 1, srcgrid%nvertlevs

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Loop through vertical coordinate

       do k = 1, dstgrid%nvertlevs

          ! Define local variables
             
          mpi_node_source      = mpi_masternode
          mpi_node_destination = k
       
          ! Send all variables to appropriate compute (slave) node
          ! (task)

          call mpi_send(dst_vort(1:nlon,1:nlat,k),(nlon*nlat),mpi_real,    &
               & mpi_node_destination,mpi_node_source,mpi_comm_world,      &
               & mpi_ierror)

       end do ! do k = 1, dstgrid%nvertlevs

       ! Loop through vertical coordinate

       do k = 1, dstgrid%nvertlevs
          
          ! Define local variables
          
          mpi_node_source      = mpi_masternode
          mpi_node_destination = k

          ! Receive all variables on root task
          
          call mpi_recv(dst_psi(1:nlon,1:nlat,k),(nlon*nlat),mpi_real,     &
               & mpi_node_destination,mpi_node_source,mpi_comm_world,      &
               & mpi_errorstatus,mpi_ierror)

          ! Initialize local variable

          count = 1

          ! Loop through meridional horizontal coordinate
          
          do j = 1, nlat
             
             ! Loop through zonal horizontal coordinate

             do i = 1, nlon

                ! Define local variable

                dst_variable(k,count) = dst_psi(i,j,k)

                ! Update counting variable

                count = count + 1

             end do ! do i = 1, nlon

          end do ! do j = 1, nlat

       end do ! do k = 1, dstgrid%nvertlevs

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%nvertlevs) then

       ! Receive all variables on appropriate compute (slave) node
       ! (task)

       call mpi_recv(vort,(nlon*nlat),mpi_real,mpi_masternode,0,           &
            & mpi_comm_world,mpi_errorstatus,mpi_ierror)

       ! Initialize local variable

       rhs = 0.0

       ! Define local variables
       
       xbdl = 0.0
       xbdr = 0.0
       ybdb = 0.0
       ybdt = 0.0

       ! Loop through meridional horizontal coordinate

       do j = 1, nlat
          
          ! Loop through zonal horizontal coordinate

          do i = 1, nlon

             ! Define local variable

             rhs(i,j) = vort(i,j)

          end do ! do i = 1, nlon

       end do ! do j = 1, nlat
       
       ! Compute local variable

       call hwscrt(xmin,xmax,nlon,0,xbdl,xbdr,ymin,ymax,nlat,1,ybdb,       &
            & ybdt,0.0,rhs,(nlon+1),perturb,ierror,work)

       ! Check local variable and proceed accordingly

       if(is_smooth) then

          ! Initialize local variable

          count = 1

          ! Loop through meridional horizontal coordinate

          do j = 1, nlat

             ! Loop through zonal horizontal coordinate

             do i = 1, nlon

                ! Define local variable

                dstgrid_var(count) = rhs(i,j)

                ! Update counting variable

                count = count + 1

             end do ! do i = 1, nlon

          end do ! do j = 1, nlat

          ! Compute local variable

          if(is_smooth) call interpolation_smoothanalysis(dstgrid_var)

          ! Initialize local variable

          count = 1

          ! Loop through meridional horizontal coordinate

          do j = 1, nlat

             ! Loop through zonal horizontal coordinate

             do i = 1, nlon

                ! Define local variable

                rhs(i,j) = dstgrid_var(count)

                ! Update counting variable

                count = count + 1

             end do ! do i = 1, nlon

          end do ! do j = 1, nlat

       end if ! if(is_smooth)

       ! Loop through meridional horizontal coordinate

       do j = 1, nlat
          
          ! Loop through zonal horizontal coordinate

          do i = 1, nlon

             ! Define local variable

             psi(i,j) = rhs(i,j)

          end do ! do i = 1, nlon

       end do ! do j = 1, nlat

       ! Send variable to master (root) node (task)

       call mpi_send(psi,(nlon*nlat),mpi_real,mpi_masternode,0,            &
            & mpi_comm_world,mpi_ierror)

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.       &
           ! dstgrid%nvertlevs)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------
    
    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dst_variable,(dstgrid%ncoords*dstgrid%nvertlevs),       &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_streamfunction

  !=======================================================================

  ! calculate_velocitypotential.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_velocitypotential(srcgrid,dstgrid,dst_variable)

    ! Define variables passed to routine

    type(interpgrid)                                                                    :: srcgrid
    type(interpgrid)                                                                    :: dstgrid
    type(grid_interface)                                                                :: grid

    ! Define variable returned by routine

    real(r_kind),              dimension(dstgrid%nvertlevs,dstgrid%ncoords)             :: dst_variable

    ! Define variables computed within routine

    real(r_double),            dimension(srcgrid%nvertlevs*srcgrid%ncoords)             :: workgrid
    real(r_kind),              dimension(nlon,nlat,dstgrid%nvertlevs)                   :: dst_divg
    real(r_kind),              dimension(nlon,nlat,dstgrid%nvertlevs)                   :: dst_chi
    real(r_kind),              dimension(nlon+1,nlat+1)                                 :: rhs
    real(r_kind),              dimension(nlon,nlat)                                     :: divg
    real(r_kind),              dimension(nlon,nlat)                                     :: chi
    real(r_kind),              dimension(nlat+1)                                        :: xbdl 
    real(r_kind),              dimension(nlat+1)                                        :: xbdr 
    real(r_kind),              dimension(nlon+1)                                        :: ybdt 
    real(r_kind),              dimension(nlon+1)                                        :: ybdb 
    real(r_kind),              dimension(10000000)                                      :: work
    real(r_kind)                                                                        :: perturb
    real(r_kind)                                                                        :: xmin
    real(r_kind)                                                                        :: xmax
    real(r_kind)                                                                        :: ymin
    real(r_kind)                                                                        :: ymax
    integer                                                                             :: ierror

    ! Define counting variables

    integer                                                                             :: i, j, k, l
    integer                                                                             :: count
    integer                                                                             :: dst_count

    !=====================================================================

    ! Initialize local variable

    dst_variable = 1.e30

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize local variables
       
       call init_constants_derived()
       call init_constants(.true.)
       perturb = 0.0

       ! Compute local variable

       call calculate_divergence(workgrid)

       ! Compute local variables

       xmin = 1.0
       xmax = dx*cos((sum(dstgrid%xlat)/real(dstgrid%ncoords)))*          &
            & ((2.0*pi*rearth_equator)/360.0)
       ymin = 1.0
       ymax = dy*cos((sum(dstgrid%xlat)/real(dstgrid%ncoords)))*111321.0

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(workgrid,(mpas_ncoords*mpas_nlevels),mpi_real,          &
         & mpi_masternode,mpi_comm_world,mpi_ierror)          
    call mpi_bcast(xmin,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(xmax,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(ymin,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(ymax,1,mpi_real,mpi_masternode,mpi_comm_world,          &
         & mpi_ierror)
    call mpi_bcast(perturb,1,mpi_real,mpi_masternode,mpi_comm_world,       &
         & mpi_ierror)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, srcgrid%nvertlevs

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then
       
          ! Loop through horizontal coordinate
                   
          do i = 1, srcgrid%ncoords
                      
             ! Define local variable
             
             srcgrid_var(i) = real(workgrid(count))
             
             ! Update counting variable
             
             count = count + 1
             
          end do ! do i = 1, srcgrid%ncoords

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)

       ! Compute local variables

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)

       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)

       ! If on master (root) node (task), define problem and broadcast
       ! variables to each slave (compute) node (task)

       if(mpi_procid .eq. mpi_masternode) then

          ! Initialize local variable

          dst_count = 1

          ! Loop through meridional horizontal coordinate

          do j = 1, nlat
       
             ! Loop through zonal horizontal coordinate
                   
             do i = 1, nlon
                      
                ! Define local variable
             
                dst_divg(i,j,k) = dstgrid_var(dst_count)
             
                ! Update counting variable
             
                dst_count = dst_count + 1
             
             end do ! do i = 1, nlon

          end do !  do j = 1, nlat

       end if ! if(mpi_procid .eq. mpi_masternode)

       ! Enable the root task to catch up from I/O and calculations

       call mpi_barrier(mpi_comm_world,mpi_ierror)

    end do ! do k = 1, srcgrid%nvertlevs

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Loop through vertical coordinate

       do k = 1, dstgrid%nvertlevs

          ! Define local variables
             
          mpi_node_source      = mpi_masternode
          mpi_node_destination = k
       
          ! Send all variables to appropriate compute (slave) node
          ! (task)

          call mpi_send(dst_divg(1:nlon,1:nlat,k),(nlon*nlat),mpi_real,    &
               & mpi_node_destination,mpi_node_source,mpi_comm_world,      &
               & mpi_ierror)

       end do ! do k = 1, dstgrid%nvertlevs

       ! Loop through vertical coordinate

       do k = 1, dstgrid%nvertlevs
          
          ! Define local variables
          
          mpi_node_source      = mpi_masternode
          mpi_node_destination = k

          ! Receive all variables on root task
          
          call mpi_recv(dst_chi(1:nlon,1:nlat,k),(nlon*nlat),mpi_real,     &
               & mpi_node_destination,mpi_node_source,mpi_comm_world,      &
               & mpi_errorstatus,mpi_ierror)

          ! Initialize local variable

          count = 1

          ! Loop through meridional horizontal coordinate
          
          do j = 1, nlat
             
             ! Loop through zonal horizontal coordinate

             do i = 1, nlon

                ! Define local variable

                dst_variable(k,count) = dst_chi(i,j,k)

                ! Update counting variable

                count = count + 1

             end do ! do i = 1, nlon

          end do ! do j = 1, nlat

       end do ! do k = 1, dstgrid%nvertlevs

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%nvertlevs) then

       ! Receive all variables on appropriate compute (slave) node
       ! (task)

       call mpi_recv(divg,(nlon*nlat),mpi_real,mpi_masternode,0,           &
            & mpi_comm_world,mpi_errorstatus,mpi_ierror)

       ! Initialize local variable

       rhs = 0.0

       ! Define local variables
       
       xbdl = 0.0
       xbdr = 0.0
       ybdb = 0.0
       ybdt = 0.0

       ! Loop through meridional horizontal coordinate

       do j = 1, nlat
          
          ! Loop through zonal horizontal coordinate

          do i = 1, nlon

             ! Define local variable

             rhs(i,j) = divg(i,j)

          end do ! do i = 1, nlon

       end do ! do j = 1, nlat
       
       ! Compute local variable

       call hwscrt(xmin,xmax,nlon,0,xbdl,xbdr,ymin,ymax,nlat,1,ybdb,       &
            & ybdt,0.0,rhs,(nlon+1),perturb,ierror,work)

       ! Check local variable and proceed accordingly

       if(is_smooth) then

          ! Initialize local variable

          count = 1

          ! Loop through meridional horizontal coordinate

          do j = 1, nlat

             ! Loop through zonal horizontal coordinate

             do i = 1, nlon

                ! Define local variable

                dstgrid_var(count) = rhs(i,j)

                ! Update counting variable

                count = count + 1

             end do ! do i = 1, nlon

          end do ! do j = 1, nlat

          ! Compute local variable

          call interpolation_smoothanalysis(dstgrid_var)

          ! Initialize local variable

          count = 1

          ! Loop through meridional horizontal coordinate

          do j = 1, nlat

             ! Loop through zonal horizontal coordinate

             do i = 1, nlon

                ! Define local variable

                rhs(i,j) = dstgrid_var(count)

                ! Update counting variable

                count = count + 1

             end do ! do i = 1, nlon

          end do ! do j = 1, nlat

       end if ! if(is_smooth)

       ! Loop through meridional horizontal coordinate

       do j = 1, nlat
          
          ! Loop through zonal horizontal coordinate

          do i = 1, nlon

             ! Define local variable

             chi(i,j) = rhs(i,j)

          end do ! do i = 1, nlon

       end do ! do j = 1, nlat

       ! Send variable to master (root) node (task)

       call mpi_send(chi,(nlon*nlat),mpi_real,mpi_masternode,0,            &
            & mpi_comm_world,mpi_ierror)

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.       &
           ! dstgrid%nvertlevs)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------
    
    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dst_variable,(dstgrid%ncoords*dstgrid%nvertlevs),       &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_velocitypotential

  !=======================================================================

  ! calculate_potentialvorticity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_potentialvorticity(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_nlevels*mpas_ncoords) :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(mpas_nlevels,mpas_ncoords) :: pressure
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords) :: theta
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords) :: dtheta_dpressure
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords) :: workgrid
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords) :: src_absvort
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords) :: src_dtheta_dpressure
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords) :: src_vort
    real(r_double),            dimension(mpas_ncoords)              :: src_cori

    ! Define counting variables

    integer                                                           :: i, j, k
    integer                                                           :: count

    !=====================================================================

    ! Initialize local variables
    
    call init_constants_derived()
    call init_constants(.true.)
    
    !---------------------------------------------------------------------    

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,    &
         & ncid=ncfileid)
    
    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Define local variable
    
    pressure = workgrid
    
    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Define local variable
    
    pressure = pressure + workgrid
    
    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'theta',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,theta)
    
    ! Close external file
    
    ncstatus = nf90_close(ncfileid)
    
    !---------------------------------------------------------------------    

    ! Compute local variables
    
    call calculate_vorticity(src_vort)
    call calculate_coriolis(src_cori)
    
    !---------------------------------------------------------------------    

    ! Initialize local variable
    
    count = 1
    
    ! Loop through vertical coordinate
    
    do k = 1, mpas_nlevels
       
       ! Loop through horizontal coordinate
       
       do i = 1, mpas_ncoords
          
          ! Check local variable and proceed accordingly
          
          if(k .eq. 1) then
             
             ! Compute local variable
             
             dtheta_dpressure(k,i) = (theta(k+1,i) - theta(k,i))/          &   
                  & (pressure(k+1,i) - pressure(k,i))
             
          end if ! if(k .eq. 1)
          
          ! Check local variable and proceed accordingly
          
          if(k .ne. 1 .and. k .ne. mpas_nlevels) then

             ! Compute local variable

             dtheta_dpressure(k,i) = (theta(k+1,i) - theta(k-1,i))/        &
                  & (pressure(k+1,i) - pressure(k-1,i))
             
          end if ! if(k .ne. 1 .and. k .ne. mpas_nlevels)

          ! Check local variable and proceed accordingly

          if(k .eq. mpas_nlevels) then

             ! Compute local variable

             dtheta_dpressure(k,i) = (theta(k,i) - theta(k-1,i))/          &
                  & (pressure(k,i) - pressure(k-1,i))

          end if ! if(k .eq. mpas_nlevels)            

          ! Define local variables

          src_absvort(count)          = src_vort(count) + src_cori(i)
          src_dtheta_dpressure(count) = dtheta_dpressure(k,i)

          ! Compute local variable

          src_variable(count) = -grav*src_absvort(count)*                  &
               & src_dtheta_dpressure(count)*1e5
             
          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_potentialvorticity

  !=======================================================================

  ! calculate_coriolis.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_coriolis(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_ncoords)                :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(mpas_ncoords)                :: src_lat

    !=====================================================================

    ! Initialize local variables
    
    call init_constants_derived()
    call init_constants(.true.)
    
    !---------------------------------------------------------------------   

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,    &
         & ncid=ncfileid)

    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'latCell',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,src_lat)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    ! Compute local variable

    src_variable = 2.0*earth_omega*(src_lat*deg2rad)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_coriolis

  !=======================================================================

  ! calculate_reflectivity_dbz.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_reflectivity_dbz(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(mpas_nlevels,mpas_ncoords)   :: workgrid_mixratio
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords)   :: workgrid_rain
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords)   :: workgrid_snow
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords)   :: workgrid_graupel
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords)   :: workgrid
    real(r_double),            dimension(mpas_nlevels,mpas_ncoords)   :: pressure
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_pressure
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_temperature
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_mixratio
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_rain
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_snow
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_graupel
    real(r_double)                                                    :: src_density
    real(r_double)                                                    :: gamma_seven
    real(r_double)                                                    :: r1
    real(r_double)                                                    :: ron
    real(r_double)                                                    :: ronv
    real(r_double)                                                    :: ron2
    real(r_double)                                                    :: son
    real(r_double)                                                    :: sonv
    real(r_double)                                                    :: gon
    real(r_double)                                                    :: gonv
    real(r_double)                                                    :: ron_min
    real(r_double)                                                    :: ron_qr0
    real(r_double)                                                    :: ron_delqr0
    real(r_double)                                                    :: ron_const1r
    real(r_double)                                                    :: ron_const2r
    real(r_double)                                                    :: rho_r
    real(r_double)                                                    :: rho_s
    real(r_double)                                                    :: rho_g
    real(r_double)                                                    :: alpha
    real(r_double)                                                    :: factor_r
    real(r_double)                                                    :: factor_s
    real(r_double)                                                    :: factor_g
    real(r_double)                                                    :: factorb_s
    real(r_double)                                                    :: factorb_g
    real(r_double)                                                    :: temp_c

    ! Define counting variables

    integer                                                           :: i, j, k
    integer                                                           :: count

    !=====================================================================

    ! Initialize local variables
    
    call init_constants_derived()
    call init_constants(.true.)

    !--------------------------------------------------------------------- 

    ! Define local variables (these follow from Thompson et al., 2004)

    r1          = dble(1.0e-15)
    ron         = dble(8.0e6)
    ron2        = dble(1.0e10)
    son         = dble(2.0e7)
    gon         = dble(5.0e7)
    ron_min     = dble(8.0e6)
    ron_qr0     = dble(0.00010)
    ron_delqr0  = dble(0.25*ron_qr0)
    ron_const1r = dble((ron2 - ron_min)*0.5)
    ron_const2r = dble((ron2 + ron_min)*0.5)
    gamma_seven = dble(720.0)
    rho_r       = dble(1000.0)
    rho_s       = dble(100.0)
    rho_g       = dble(400.0)
    alpha       = dble(0.224)

    !---------------------------------------------------------------------    

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,    &
         & ncid=ncfileid)
    
    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Define local variable
    
    pressure = workgrid
    
    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    
    ! Define local variable
    
    pressure = pressure + workgrid

    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'qv',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid_mixratio)

    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'qr',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid_rain)
    
    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'qs',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid_snow)

    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'qg',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid_graupel)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------   

    ! Compute local variable

    call calculate_temperature(src_temperature)

    ! Rescale local variables accordingly

    where(src_temperature .lt. 273.15) src_snow = src_rain
    where(src_temperature .lt. 273.15) src_rain = dble(0.0)

    !--------------------------------------------------------------------- 

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, mpas_nlevels

       ! Loop through horizontal coordinate

       do i = 1, mpas_ncoords

          ! Define local variables

          src_pressure(count) = pressure(k,i)
          src_mixratio(count) = workgrid_mixratio(k,i)
          src_rain(count)     = workgrid_rain(k,i)
          src_snow(count)     = workgrid_snow(k,i)
          src_graupel(count)  = workgrid_graupel(k,i)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    ! Update local variables accordingly

    src_mixratio = dble(max(real(src_mixratio),0.0))
    src_rain     = dble(max(real(src_rain),0.0))
    src_snow     = dble(max(real(src_snow),0.0))
    src_graupel  = dble(max(real(src_graupel),0.0))

    ! Compute local variables

    factor_r = dble(gamma_seven*dble(1.0e18)*(1.0/(dble(pi)*rho_r))**1.75)
    factor_s = dble(gamma_seven*dble(1.0e18)*(1.0/(dble(pi)*rho_s))**1.75  &
         & *(rho_s/rho_r)**2.0*alpha)
    factor_g = dble(gamma_seven*dble(1.0e18)*(1.0/(dble(pi)*rho_s))**1.75  &
         & *(rho_g/rho_r)**2.0*alpha)

    !---------------------------------------------------------------------   

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, mpas_nlevels

       ! Loop through horizontal coordinate

       do i = 1, mpas_ncoords

          ! Compute local variable

          src_density = src_pressure(count)/(287.04*                       &
               & src_temperature(count)*(1.0 +                             &
               & 0.608*src_mixratio(count)))

          ! Check local variable and proceed accordingly

          if(src_temperature(count) .gt. 273.15) then

             ! Compute local variables

             factorb_s = factor_s/alpha
             factorb_g = factor_g/alpha

          else   ! if(src_temperature(count) .gt. 273.15)

             ! Define local variables

             factorb_s = factor_s
             factorb_g = factor_g

          end if ! if(src_temperature(count) .gt. 273.15)

          ! Define local variables

          temp_c = dble(amin1(-0.001,real(src_temperature(count))-275.15))
          sonv   = dble(amin1(2.0e8,2.0e6*exp(-0.12*real(temp_c))))
          gonv   = gon
          ronv   = ron2

          ! Check local variable and proceed accordingly

          if(src_graupel(count) .gt. r1) then
             
             ! Compute local variables

             gonv = dble(2.38*(pi*rho_g/(src_density*src_graupel(count)))  &
                  & **0.92)
             gonv = dble(max(1.0e4,min(real(gonv),real(gon))))

          end if ! if(src_graupel(count) .gt. r1)

          ! Check local variable and proceed accordingly

          if(src_rain(count) .gt. r1) then
          
             ! Compute local variables

             ronv = dble(ron_const1r*tanh((real(ron_qr0)-                  &
                  & real(src_rain(count)))/ron_delqr0) + ron_const2r)

          end if ! if(src_rain(count) .gt. r1)

          ! Compute local variable

          src_variable(count) = dble(factor_r*(src_density*                &
               & src_rain(count))**1.75/ronv**0.75 +                       &
               & factorb_s*(src_density*src_snow(count))**1.75/sonv**0.75  &
               & + factorb_g*(src_density*src_graupel(count))**1.75/       &
               & gonv**0.75)

          ! Rescale local variable accordingly

          src_variable(count) = dble(max(real(src_variable(count)),0.001))
          src_variable(count) = dble(10.0*log10(                           &
               & real(src_variable(count))))
 
          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_reflectivity_dbz

  !=======================================================================

  ! calculate_reflectivity_composite.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_reflectivity_composite(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_ncoords)                :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(mpas_nlevels,mpas_ncoords)   :: workgrid
    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_dbz

    ! Define counting variables

    integer                                                           :: i, j, k
    integer                                                           :: count

    !=====================================================================

    ! Initialize local variables
    
    call init_constants_derived()
    call init_constants(.true.)

    !---------------------------------------------------------------------  

    ! Compute local variable

    call calculate_reflectivity_dbz(src_dbz)

    !--------------------------------------------------------------------- 

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, mpas_nlevels

       ! Loop through horizontal coordinate

       do i = 1, mpas_ncoords

          ! Define local variable

          workgrid(k,i) = src_dbz(count)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, mpas_nlevels

    ! Loop through horizontal coordinate

    do i = 1, mpas_ncoords

       ! Define local variable

       src_variable(i) = maxval(workgrid(1:mpas_nlevels,i))

    end do ! do i = 1, mpas_ncoords

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_reflectivity_composite

  !=======================================================================

  ! calculate_height.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_height(src_variable)

    ! Define variable returned by routine

    real(r_double),            dimension(mpas_nlevels*mpas_ncoords)   :: src_variable

    ! Define variables computed within routine

    real(r_double),            dimension(mpas_nlevels+1,mpas_ncoords) :: src_height

    ! Define counting variables

    integer                                                           :: i, j, k 
    integer                                                           :: count

    !=====================================================================

    ! Initialize local variables
    
    call init_constants_derived()
    call init_constants(.true.)

    ! Open external file
    
    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,    &
         & ncid=ncfileid)

    ! Define and ingest local variables
    
    ncstatus = nf90_inq_varid(ncfileid,'zgrid',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,src_height)

    ! Close external file
    
    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1

    ! Loop through vertical coordinate

    do k = 1, mpas_nlevels

       ! Loop through horizontal coordinate

       do i = 1, mpas_ncoords

          ! Compute local variable

          src_variable(count) = (src_height(k,i) + src_height((k+1),i))     &
               & /2.0

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, mpas_ncoords

    end do ! do k = 1, (mpas_nlevels + 1)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine calculate_height

  !=======================================================================

  ! calculate_clipvariable.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_clipvariable(xdim,zdim,grid)

    ! Define array dimension variables

    integer                                                  :: xdim
    integer                                                  :: zdim

    ! Define variables passed to routine

    real(r_kind),            dimension(xdim*zdim)            :: grid

    ! Define variables computed within routine

    real(r_single)                                           :: clipval

    !=====================================================================

    ! Define local variable

    clipval = tiny(grid(1))

    ! Update local variable accordingly

    where(grid .lt. clipval) grid = clipval

    !=====================================================================

    ! Return calculated variables

    return

    !=====================================================================

  end subroutine calculate_clipvariable

  !=======================================================================

end module variable_interface
