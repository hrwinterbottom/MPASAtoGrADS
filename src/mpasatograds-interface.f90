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

module mpastograds_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use grads_interface
  use interpolation_interface
  use mpi_interface
  use namelist
  use variable_interface

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

  ! Define global variables

  real(r_kind),                    dimension(:,:,:), allocatable :: mpas_pressure
  real(r_kind),                    dimension(:,:,:), allocatable :: mpas_height
  real(r_kind),                    dimension(:,:,:), allocatable :: mpas_psfc

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! mpastograds.f90:

  !-----------------------------------------------------------------------

  subroutine mpastograds()

    ! Define variables computed within routine

    type(variable_info),      dimension(:),      allocatable :: mpas
    type(grid_interface)                                     :: grid
    type(interpgrid)                                         :: srcgrid
    type(interpgrid)                                         :: dstgrid

    !=====================================================================

    ! Ingest external mpastograds.input

    call namelistparams()

    !---------------------------------------------------------------------

    ! Define local variables

    call interpolation_initialize_define_grid()
    call mpas_dimensions()

    ! Compute local variables

    call mpas_interpolation_grids_initialize(grid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%atmos_mass_ncells,1,mpi_integer,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%gauss_mass_ncells,1,mpi_integer,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(barnes_npasses,1,mpi_integer,mpi_masternode,            &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(barnes_nneighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%atmos_mass_lon,grid%atmos_mass_ncells,mpi_double,  &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%atmos_mass_lat,grid%atmos_mass_ncells,mpi_double,  &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%gauss_mass_lon,grid%gauss_mass_ncells,mpi_real,    &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%gauss_mass_lat,grid%gauss_mass_ncells,mpi_real,    &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(barnes_distance_threshold,1,mpi_real,mpi_masternode,    & 
         & mpi_comm_world,mpi_ierror)

    ! Define local variables

    srcgrid%ncoords   = grid%atmos_mass_ncells
    srcgrid%nvertlevs = mpas_nlevels
    srcgrid%npasses   = barnes_npasses
    srcgrid%neighbors = barnes_nneighbors
    dstgrid%ncoords   = grid%gauss_mass_ncells
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

    srcgrid%xlong      = real(grid%atmos_mass_lon)
    srcgrid%xlat       = real(grid%atmos_mass_lat)
    srcgrid%distthresh = barnes_distance_threshold
    dstgrid%xlong      = grid%gauss_mass_lon
    dstgrid%xlat       = grid%gauss_mass_lat
    dstgrid%distthresh = barnes_distance_threshold

    ! Compute local variables

    call interpolation_define_kdtree_mpi(srcgrid)
    call interpolation_define_kdtree_mpi(dstgrid)
    call define_scaling_coefficients(srcgrid)
    call define_scaling_coefficients(dstgrid)
    call interpolation_define_weights_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*dstgrid%neighbors*     &
         & dstgrid%npasses),mpi_real,mpi_masternode,mpi_comm_world,        &
         & mpi_ierror)
    call mpi_bcast(dstgrid%grdnbors,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_integer,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Define local variables

    grid%atmos_mass_2_gauss_mass_weights   = dstgrid%weights
    grid%atmos_mass_2_gauss_mass_neighbors = dstgrid%grdnbors

    ! Deallocate memory for local variables

    call interpolation_cleanup_grid(srcgrid)
    call interpolation_cleanup_grid(dstgrid)
    call interpolation_cleanup_task_balance(srcgrid)
    call interpolation_cleanup_task_balance(dstgrid)

    ! Define local variables

    srcgrid%ncoords   = grid%gauss_mass_ncells
    srcgrid%npasses   = barnes_npasses
    srcgrid%neighbors = barnes_nneighbors
    dstgrid%ncoords   = grid%atmos_mass_ncells
    dstgrid%npasses   = barnes_npasses
    dstgrid%neighbors = barnes_nneighbors

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(srcgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%npasses,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%ncoords,1,mpi_integer,mpi_masternode,           &
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

    srcgrid%xlong      = grid%gauss_mass_lon
    srcgrid%xlat       = grid%gauss_mass_lat
    srcgrid%distthresh = barnes_distance_threshold
    dstgrid%xlong      = real(grid%atmos_mass_lon)
    dstgrid%xlat       = real(grid%atmos_mass_lat)
    dstgrid%distthresh = barnes_distance_threshold

    ! Compute local variables

    call interpolation_define_kdtree_mpi(srcgrid)
    call interpolation_define_kdtree_mpi(dstgrid)
    call define_scaling_coefficients(srcgrid)
    call define_scaling_coefficients(dstgrid)
    call interpolation_define_weights_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*dstgrid%neighbors*     &
         & dstgrid%npasses),mpi_real,mpi_masternode,mpi_comm_world,        &
         & mpi_ierror)
    call mpi_bcast(dstgrid%grdnbors,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_integer,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Define local variables

    grid%gauss_mass_2_atmos_mass_weights   = dstgrid%weights
    grid%gauss_mass_2_atmos_mass_neighbors = dstgrid%grdnbors

    ! Deallocate memory for local variables

    call interpolation_cleanup_grid(srcgrid)
    call interpolation_cleanup_grid(dstgrid)
    call interpolation_cleanup_task_balance(srcgrid)
    call interpolation_cleanup_task_balance(dstgrid)

    !---------------------------------------------------------------------

    ! Define local variable
       
    call variable_interface_num_var()
    
    ! Initialize local variables
    
    call variable_interface_initialize(mpas)
    var_counter = 0
    
    ! Define local variables
    
    call variable_interface_isolevel(mpas)
    call variable_interface_profile(mpas)
       
    ! Initialize local variables
       
    call variable_interface_derived_initialize(mpas)

    ! Compute local variables

    call variable_interface_process(mpas,srcgrid,dstgrid,grid)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Write external files

       call grads_interface_descriptor(mpas)
       call grads_interface_binary(mpas)

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Deallocate memory for local variable

    call variable_interface_cleanup(mpas)

    !=====================================================================

  end subroutine mpastograds

  !=======================================================================

  subroutine mpas_dimensions()

    !=====================================================================

    ! Define local variables

    ncstatus               = nf90_open(path=trim(variable_filename),       &
         & mode=nf90_nowrite,ncid=ncfileid) 
    ncstatus = nf90_inq_dimid(ncfileid,'nCells',ncdimid)
    ncstatus               = nf90_inquire_dimension(ncfileid,ncdimid,      &
         & len=mpas_ncoords)
    ncstatus = nf90_inq_dimid(ncfileid,'nVertLevels',ncdimid)
    ncstatus               = nf90_inquire_dimension(ncfileid,ncdimid,      &
         & len=mpas_nlevels)
    ncstatus               = nf90_close(ncfileid)

    !=====================================================================

  end subroutine mpas_dimensions

  !=======================================================================

  ! mpas_interpolation_grids_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine mpas_interpolation_grids_initialize(grid)

    ! Define variables computed within routine

    type(grid_interface)                                                 :: grid

    ! Define counting variables

    integer                                                              :: i, j, k
    integer                                                              :: count

    !=====================================================================

    ! Define local variables

    ncstatus               = nf90_open(path=trim(variable_filename),       &
         & mode=nf90_nowrite,ncid=ncfileid) 
    ncstatus = nf90_inq_dimid(ncfileid,'nCells',ncdimid)
    ncstatus               = nf90_inquire_dimension(ncfileid,ncdimid,      &
         & len=grid%atmos_mass_ncells)
    ncstatus               = nf90_close(ncfileid)
    grid%gauss_mass_ncells = nlon*nlat

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(grid%atmos_mass_lon))                               &
         & allocate(grid%atmos_mass_lon(grid%atmos_mass_ncells))
    if(.not. allocated(grid%atmos_mass_lat))                               &
         & allocate(grid%atmos_mass_lat(grid%atmos_mass_ncells))
    if(.not. allocated(grid%gauss_mass_lon))                               &
         & allocate(grid%gauss_mass_lon(grid%gauss_mass_ncells))
    if(.not. allocated(grid%gauss_mass_lat))                               &
         & allocate(grid%gauss_mass_lat(grid%gauss_mass_ncells))

    !---------------------------------------------------------------------

    ! Define local variables

    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)
    ncstatus = nf90_inq_varid(ncfileid,'lonCell',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%atmos_mass_lon)
    ncstatus = nf90_inq_varid(ncfileid,'latCell',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,grid%atmos_mass_lat)
    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(grid%atmos_mass_2_gauss_mass_weights))              &
         & allocate(grid%atmos_mass_2_gauss_mass_weights(                  &
         & grid%gauss_mass_ncells,barnes_nneighbors,barnes_npasses))
    if(.not. allocated(grid%gauss_mass_2_atmos_mass_weights))              &
         & allocate(grid%gauss_mass_2_atmos_mass_weights(                  &
         & grid%atmos_mass_ncells,barnes_nneighbors,barnes_npasses))
    if(.not. allocated(grid%atmos_mass_2_gauss_mass_neighbors))            &
         & allocate(grid%atmos_mass_2_gauss_mass_neighbors(                &
         & grid%gauss_mass_ncells,barnes_nneighbors))
    if(.not. allocated(grid%gauss_mass_2_atmos_mass_neighbors))            &
         & allocate(grid%gauss_mass_2_atmos_mass_neighbors(                &
         & grid%atmos_mass_ncells,barnes_nneighbors))

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(grid_xlong)) allocate(grid_xlong(grid_xdim))
    if(.not. allocated(grid_xlat))  allocate(grid_xlat(grid_ydim))

    ! Allocate memory for local variable
    
    if(.not. allocated(grid%gauss_mass_lon))                               &
         & allocate(grid%gauss_mass_lon(grid%gauss_mass_ncells)) 
    if(.not. allocated(grid%gauss_mass_lat))                               &
         & allocate(grid%gauss_mass_lat(grid%gauss_mass_ncells))
    
    ! Compute local variables

    call interpolation_define_grid()

    ! Initialize counting variable

    count = 1

    ! Loop through meridional horizontal coordinate

    do j = 1, nlat

       ! Loop through zonal horizontal coordinate

       do i = 1, nlon

          ! Define local variables

          grid%gauss_mass_lon(count) = grid_xlong(i)
          grid%gauss_mass_lat(count) = grid_xlat(j)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, nlon
 
    end do ! do j = 1, nlat

    ! Deallocate memory for local variables

    if(allocated(grid_xlong)) deallocate(grid_xlong)
    if(allocated(grid_xlat))  deallocate(grid_xlat)

    !=====================================================================

  end subroutine mpas_interpolation_grids_initialize

  !=======================================================================

  ! mpas_interpolation_grid_cleanup.f90:

  !-----------------------------------------------------------------------

  subroutine mpas_interpolation_grids_cleanup(grid)

    ! Define variables passed to routine

    type(grid_interface)                                                 :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%atmos_mass_lon))                                     &
         & deallocate(grid%atmos_mass_lon)
    if(allocated(grid%atmos_mass_lat))                                     &
         & deallocate(grid%atmos_mass_lat)
    if(allocated(grid%gauss_mass_lon))                                     &
         & deallocate(grid%gauss_mass_lon)
    if(allocated(grid%gauss_mass_lat))                                     &
         & deallocate(grid%gauss_mass_lat)
    if(allocated(grid%atmos_mass_2_gauss_mass_weights))                    &
         & deallocate(grid%atmos_mass_2_gauss_mass_weights)
    if(allocated(grid%gauss_mass_2_atmos_mass_weights))                    &
         & deallocate(grid%gauss_mass_2_atmos_mass_weights)
    if(allocated(grid%atmos_mass_2_gauss_mass_neighbors))                  &
         & deallocate(grid%atmos_mass_2_gauss_mass_neighbors)
    if(allocated(grid%gauss_mass_2_atmos_mass_neighbors))                  &
         & deallocate(grid%gauss_mass_2_atmos_mass_neighbors)

    !=====================================================================

  end subroutine mpas_interpolation_grids_cleanup

  !=======================================================================

end module mpastograds_interface
