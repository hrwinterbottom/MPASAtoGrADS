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

!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License along
!    with this program; if not, write to the Free Software Foundation, Inc.,
!    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

module interpolation_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use kdtree2_module
  use mpi_interface
  use namelist
  use netcdf
  use sphere

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

  ! Define all data and structure types for routine; these variables
  ! are variables required by the subroutines within this module

  type interpgrid
     type(kdtree2),        pointer                               :: kdtree_grid
     real(r_kind),                 dimension(:,:,:), allocatable :: weights
     real(r_kind),                 dimension(:,:),   allocatable :: anlysvar
     real(r_kind),                 dimension(:,:),   allocatable :: pressure_profile
     real(r_kind),                 dimension(:,:),   allocatable :: height_profile
     real(r_kind),                 dimension(:,:),   allocatable :: thetae_profile
     real(r_kind),                 dimension(:),     allocatable :: surface_pressure
     real(r_kind),                 dimension(:),     allocatable :: surface_height
     real(r_kind),                 dimension(:),     allocatable :: xlong
     real(r_kind),                 dimension(:),     allocatable :: xlat
     real(r_kind),                 dimension(:),     allocatable :: mask
     real(r_kind),                 dimension(:),     allocatable :: scutoff
     real(r_kind)                                                :: distthresh
     real(r_kind)                                                :: dx
     real(r_kind)                                                :: dy
     real(r_kind),                 dimension(:,:),   allocatable :: grdloc
     integer,                      dimension(:,:),   allocatable :: grdnbors
     integer,                      dimension(:),     allocatable :: mpi_count_begin
     integer,                      dimension(:),     allocatable :: mpi_count_end
     integer                                                     :: mpi_maxprocid
     integer                                                     :: mtrunc
     integer                                                     :: ncoords
     integer                                                     :: nvertlevs
     integer                                                     :: npasses
     integer                                                     :: neighbors
  end type interpgrid

  type interppres
     real(r_kind),                 dimension(:,:),   allocatable :: profile_pres
     real(r_kind),                 dimension(:,:),   allocatable :: profile_var
     integer                                                     :: ncoords
     integer                                                     :: nsig
  end type interppres

  type interphght
     real(r_kind),                 dimension(:,:),   allocatable :: profile_hght
     real(r_kind),                 dimension(:,:),   allocatable :: profile_var
     integer                                                     :: ncoords
     integer                                                     :: nsig
  end type interphght

  type interpthetae
     real(r_kind),                 dimension(:,:),   allocatable :: profile_thetae
     real(r_kind),                 dimension(:,:),   allocatable :: profile_var
     integer                                                     :: ncoords
     integer                                                     :: nsig
  end type interpthetae

  type grid_interface
     real(r_double),               dimension(:),     allocatable :: atmos_mass_lon
     real(r_double),               dimension(:),     allocatable :: atmos_mass_lat
     real(r_kind),                 dimension(:),     allocatable :: gauss_mass_lon
     real(r_kind),                 dimension(:),     allocatable :: gauss_mass_lat
     real(r_kind),                 dimension(:,:,:), allocatable :: atmos_mass_2_gauss_mass_weights
     real(r_kind),                 dimension(:,:,:), allocatable :: gauss_mass_2_atmos_mass_weights
     integer,                      dimension(:,:),   allocatable :: atmos_mass_2_gauss_mass_neighbors
     integer,                      dimension(:,:),   allocatable :: gauss_mass_2_atmos_mass_neighbors 
     integer                                                     :: atmos_mass_ncells
     integer                                                     :: gauss_mass_ncells
  end type grid_interface

  ! Define global variables

  real(r_kind),                   dimension(:),      allocatable :: grid_xlong
  real(r_kind),                   dimension(:),      allocatable :: grid_xlat
  real(r_kind),                   dimension(:),      allocatable :: srcgrid_var
  real(r_kind),                   dimension(:),      allocatable :: dstgrid_var
  real(r_kind)                                                   :: dx
  real(r_kind)                                                   :: dy
  real(r_kind)                                                   :: rlon_min
  real(r_kind)                                                   :: rlon_max
  real(r_kind)                                                   :: rlat_min
  real(r_kind)                                                   :: rlat_max
  integer                                                        :: nlon
  integer                                                        :: nlat
  integer                                                        :: grid_xdim
  integer                                                        :: grid_ydim
  integer                                                        :: grid_nlevs

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! interpolation_initialize_define_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_define_grid()

    !=====================================================================

    ! Define local variables accordingly

    if(grid_mtrunc .eq. 62) then
       
       ! Define local variables
          
       nlon     = 192
       nlat     = 94
       dx       = 1.875000
       dy       = 1.875000
       rlon_min = 0.000
       rlon_max = 358.125
       rlat_min = -88.542
       rlat_max = 88.542
       
    end if ! if(grid_mtrunc .eq. 62)

    ! Define local variables accordingly

    if(grid_mtrunc .eq. 126) then
          
       ! Define local variables
          
       nlon     = 384
       nlat     = 190
       dx       = 0.937500
       dy       = 0.937500
       rlon_min = 0.000
       rlon_max = 359.0625
       rlat_min = -89.277
       rlat_max = 89.277
          
    end if ! if(grid_mtrunc .eq. 126)
          
    ! Define local variables accordingly

    if(grid_mtrunc .eq. 254)  then
          
       ! Define local variables
          
       nlon     = 768
       nlat     = 384
       dx       = 0.469000
       dy       = 0.469000
       rlon_min = 0.000
       rlon_max = 359.531
       rlat_min = -89.642
       rlat_max = 89.642
       
    end if ! if(grid_mtrunc .eq. 254)
   
    ! Define local variables accordingly

    if(grid_mtrunc .eq. 382)  then
       
       ! Define local variables
       
       nlon     = 1152
       nlat     = 576
       dx       = 0.313000
       dy       = 0.313000
       rlon_min = 0.000
       rlon_max = 359.687
       rlat_min = -89.761
       rlat_max = 89.761
       
    end if ! if(grid_mtrunc .eq. 382)
   
    ! Define local variables accordingly
    
    if(grid_mtrunc .eq. 574)  then
       
       ! Define local variables
       
       nlon     = 1760
       nlat     = 880
       dx       = 0.205000
       dy       = 0.205000
       rlon_min = 0.000
       rlon_max = 359.795
       rlat_min = -89.844
       rlat_max = 89.844
          
    end if ! if(grid_mtrunc .eq. 574)

    ! Define local variables

    grid_xdim = nlon
    grid_ydim = nlat

    !---------------------------------------------------------------------

    ! Print message to user

    if(mpi_procid .eq. mpi_masternode .and. debug) write(6,500) nlon, nlat
    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_INITIALIZE_DEFINE_GRID: nlon, nlat = ', i6,1x,   &
         & i6)

    !=====================================================================

  end subroutine interpolation_initialize_define_grid

  !=======================================================================

  ! interpolation_initialize_task_balance.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_task_balance(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    ! Define variables computed within routine

    integer                                                              :: mpi_count_interval

    ! Define counting variables

    integer                                                              :: i, j, k, l
    integer                                                              :: count

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%mpi_count_begin))                              &
         & allocate(grid%mpi_count_begin((mpi_nprocs)))
    if(.not. allocated(grid%mpi_count_end))                                &
         & allocate(grid%mpi_count_end((mpi_nprocs)))

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize counting variable

       count = 1

       ! Compute local variable

       mpi_count_interval = grid%ncoords/(mpi_nprocs - 1)

       ! Initialize local variables

       grid%mpi_count_begin    = 0
       grid%mpi_count_end      = 0
       grid%mpi_count_begin(1) = 1
       grid%mpi_count_end(1)   = grid%mpi_count_begin(1) +                 &
            & mpi_count_interval

       ! Loop through total number of processors

       do l = 2, mpi_nprocs - 1

          ! Define local variables

          grid%mpi_count_begin(l) = grid%mpi_count_end(l-1) + 1
          grid%mpi_count_end(l)   = grid%mpi_count_begin(l) +              &
               & mpi_count_interval

          ! Check local variables and proceed accordingly

          if(grid%mpi_count_begin(l) .gt. grid%ncoords) then

             ! Define local variables

             grid%mpi_count_begin(l) = grid%ncoords
             grid%mpi_count_end(l)   = grid%ncoords
             grid%mpi_maxprocid      = l

             ! Define exit from loop

             goto 1000

          end if ! if(grid%mpi_count_begin(l) .gt. grid%ncoords)

          ! Check local variables and proceed accordingly

          if(grid%mpi_count_end(l) .gt. grid%ncoords) then

             ! Define local variables

             grid%mpi_count_end(l) = grid%ncoords
             grid%mpi_maxprocid    = l

             ! Define exit from loop

             goto 1000

          end if ! if(grid%mpi_count_end(l) .gt. grid%ncoords)

       end do ! do l = 2, (mpi_nprocs - 1)

       ! Define exit from loop

1000   continue

       ! Loop through local variable and proceed accordingly

       do l = 1, grid%mpi_maxprocid

          ! Print message to user

          if(debug) write(6,500) grid%ncoords, l, grid%mpi_count_begin(l), &
               & grid%mpi_count_end(l)

       end do ! do l = 1, grid%mpi_maxprocid

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%mpi_count_begin,mpi_nprocs,mpi_integer,            &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%mpi_count_end,mpi_nprocs,mpi_integer,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%mpi_maxprocid,1,mpi_integer,mpi_masternode,        &
         & mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_INITIALIZE_TASK_BALANCE: (grid size/',           &
         & 'task ID/tile min/tile max) : ', i9, i6, i9, i9)

    !=====================================================================

  end subroutine interpolation_initialize_task_balance

  !=======================================================================

  ! interpolation_cleanup_task_balance.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_task_balance(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    !=====================================================================

    ! Deallocate memory for local variable

    if(allocated(grid%mpi_count_begin))                                    &
         & deallocate(grid%mpi_count_begin)
    if(allocated(grid%mpi_count_end))                                      &
         & deallocate(grid%mpi_count_end)

    !=====================================================================

  end subroutine interpolation_cleanup_task_balance

  !=======================================================================

  ! interpolation_initialize_interppres.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_interppres(grid)

    ! Define variables passed to routine

    type(interppres)                                         :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%profile_pres))                                 &
         & allocate(grid%profile_pres(grid%ncoords,grid%nsig))
    if(.not. allocated(grid%profile_var))                                  &
         & allocate(grid%profile_var(grid%ncoords,grid%nsig))

    !=====================================================================

  end subroutine interpolation_initialize_interppres

  !=======================================================================

  ! interpolation_cleanup_interppres.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_interppres(grid)

    ! Define variables passed to routine

    type(interppres)                                         :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%profile_pres)) deallocate(grid%profile_pres)
    if(allocated(grid%profile_var))  deallocate(grid%profile_var)

    !=====================================================================

  end subroutine interpolation_cleanup_interppres

  !=======================================================================

  ! interpolation_initialize_interphght.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_interphght(grid)

    ! Define variables passed to routine

    type(interphght)                                         :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%profile_hght))                                 &
         & allocate(grid%profile_hght(grid%ncoords,grid%nsig))
    if(.not. allocated(grid%profile_var))                                  &
         & allocate(grid%profile_var(grid%ncoords,grid%nsig))

    !=====================================================================

  end subroutine interpolation_initialize_interphght

  !=======================================================================

  ! interpolation_cleanup_interphght.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_interphght(grid)

    ! Define variables passed to routine

    type(interphght)                                         :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%profile_hght)) deallocate(grid%profile_hght)
    if(allocated(grid%profile_var))  deallocate(grid%profile_var)

    !=====================================================================

  end subroutine interpolation_cleanup_interphght

  !=======================================================================

  ! interpolation_initialize_interpthetae.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_interpthetae(grid)

    ! Define variables passed to routine

    type(interpthetae)                                       :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%profile_thetae))                               &
         & allocate(grid%profile_thetae(grid%ncoords,grid%nsig))
    if(.not. allocated(grid%profile_var))                                  &
         & allocate(grid%profile_var(grid%ncoords,grid%nsig))

    !=====================================================================

  end subroutine interpolation_initialize_interpthetae

  !=======================================================================

  ! interpolation_cleanup_interpthetae.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_interpthetae(grid)

    ! Define variables passed to routine

    type(interpthetae)                                       :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%profile_thetae)) deallocate(grid%profile_thetae)
    if(allocated(grid%profile_var))    deallocate(grid%profile_var)

    !=====================================================================

  end subroutine interpolation_cleanup_interpthetae

  !=======================================================================

  ! interpolation_initialize_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_grid(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%grdloc))                                             &
         & deallocate(grid%grdloc)
    if(allocated(grid%weights))                                            &
         & deallocate(grid%weights)
    if(allocated(grid%xlong))                                              &
         & deallocate(grid%xlong)
    if(allocated(grid%xlat))                                               &
         & deallocate(grid%xlat)
    if(allocated(grid%grdnbors))                                           &
         & deallocate(grid%grdnbors)
    if(allocated(grid%scutoff))                                            &
         & deallocate(grid%scutoff)
    if(allocated(grid%surface_pressure))                                   &
         & deallocate(grid%surface_pressure)
    if(allocated(grid%pressure_profile))                                   &
         & deallocate(grid%pressure_profile)
    if(allocated(grid%height_profile))                                     &
         & deallocate(grid%height_profile)
    if(allocated(grid%surface_height))                                     &
         & deallocate(grid%surface_height)
    if(allocated(grid%thetae_profile))                                     &
         & deallocate(grid%thetae_profile)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(grid%grdloc))                                       &
         & allocate(grid%grdloc(3,grid%ncoords))
    if(.not. allocated(grid%weights))                                      &
         & allocate(grid%weights(grid%ncoords,grid%neighbors,              &
         & grid%npasses))
    if(.not. allocated(grid%xlong))                                        &
         & allocate(grid%xlong(grid%ncoords))
    if(.not. allocated(grid%xlat))                                         &
         & allocate(grid%xlat(grid%ncoords))
    if(.not. allocated(grid%grdnbors))                                     &
         & allocate(grid%grdnbors(grid%ncoords,grid%neighbors))
    if(.not. allocated(grid%scutoff))                                      &
         & allocate(grid%scutoff(grid%npasses))
    if(.not. allocated(grid%surface_pressure))                             &
         & allocate(grid%surface_pressure(grid%ncoords))
    if(.not. allocated(grid%pressure_profile))                             &
         & allocate(grid%pressure_profile(grid%nvertlevs,grid%ncoords))
    if(.not. allocated(grid%height_profile))                               &
         & allocate(grid%height_profile(grid%nvertlevs,grid%ncoords))
    if(.not. allocated(grid%surface_height))                               &
         & allocate(grid%surface_height(grid%ncoords))
    if(.not. allocated(grid%thetae_profile))                               &
         & allocate(grid%thetae_profile(grid%nvertlevs,grid%ncoords))

    !=====================================================================

  end subroutine interpolation_initialize_grid

  !=======================================================================

  ! interpolation_cleanup_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_grid(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%grdloc))                                             &
         & deallocate(grid%grdloc)
    if(allocated(grid%weights))                                            &
         & deallocate(grid%weights)
    if(allocated(grid%xlong))                                              &
         & deallocate(grid%xlong)
    if(allocated(grid%xlat))                                               &
         & deallocate(grid%xlat)
    if(allocated(grid%grdnbors))                                           &
         & deallocate(grid%grdnbors)
    if(allocated(grid%scutoff))                                            &
         & deallocate(grid%scutoff)
    if(allocated(grid%surface_pressure))                                   &
         & deallocate(grid%surface_pressure)
    if(allocated(grid%pressure_profile))                                   &
         & deallocate(grid%pressure_profile)
    if(allocated(grid%height_profile))                                     &
         & deallocate(grid%height_profile)
    if(allocated(grid%surface_height))                                     &
         & deallocate(grid%surface_height)
    if(allocated(grid%thetae_profile))                                     &
         & deallocate(grid%thetae_profile)
    if(mpi_procid .eq. mpi_masternode .and. associated(grid%kdtree_grid))  &
         & call kdtree2_destroy(grid%kdtree_grid)

    !=====================================================================

  end subroutine interpolation_cleanup_grid

  !=======================================================================

  ! interpolation_define_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_grid()

    ! Define variables computed within routine

    real(r_kind),             dimension(:),      allocatable :: grid_slat
    real(r_kind),             dimension(:),      allocatable :: grid_wlat
    real(r_kind),             dimension(:),      allocatable :: workgrid

    ! Define counting variables

    integer                                                  :: i, j, k
    integer                                                  :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
    
    if(mpi_procid .eq. mpi_masternode) then

       ! Define local variables

       grid_xdim = nlon
       grid_ydim = nlat

       ! Allocate memory for local variables

       if(.not. allocated(grid_slat)) allocate(grid_slat(grid_ydim))
       if(.not. allocated(grid_wlat)) allocate(grid_wlat(grid_ydim))
       if(.not. allocated(workgrid))  allocate(workgrid(grid_ydim))

       ! Compute local variables

       call gausslat(grid_ydim,grid_slat,grid_wlat)

       ! Define local variable

       grid_xlat = acos(grid_slat) - pi/2.0
       workgrid  = grid_xlat

       ! Loop through meridional horizontal coordinate

       do j = 1, grid_ydim

          ! Define local variable

          workgrid(j) = grid_xlat(grid_ydim - j + 1)

       end do ! do j = 1, grid_ydim

       ! Define local variable

       grid_xlat = workgrid

       ! Loop through zonal horizontal coordinate

       do i = 1, grid_xdim

          ! Compute local variable

          grid_xlong(i) = (rlon_min + (i-1)*dx)*deg2rad

       end do ! do i = 1, grid_xdim

       ! Deallocate memory for local variables

       if(allocated(grid_slat)) deallocate(grid_slat)
       if(allocated(grid_wlat)) deallocate(grid_wlat)
       if(allocated(workgrid))  deallocate(workgrid)

    end if ! if(mpi_procid .eq. mpi_masternode) 

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid_xlong,nlon,mpi_real,mpi_masternode,mpi_comm_world, &
         & mpi_ierror)
    call mpi_bcast(grid_xlat,nlat,mpi_real,mpi_masternode,mpi_comm_world,  &
         & mpi_ierror)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_DEFINE_GRID: nlon, nlat = ', i6,1x,i6)

    !=====================================================================

  end subroutine interpolation_define_grid

  !=======================================================================

  ! interpolation_define_mask.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_mask(srcgrid,dstgrid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: srcgrid
    type(interpgrid)                                         :: dstgrid

    ! Define variables computed in routine

    ! Define variables computed within routine

    type(kdtree2_result),      dimension(1)                  :: sresults
    real(r_double),            dimension(srcgrid%ncoords)    :: workgrid_seaice
    real(r_double),            dimension(srcgrid%ncoords)    :: workgrid_xland
    integer                                                  :: ncfileid
    integer                                                  :: ncvarid
    integer                                                  :: ncdimid
    integer                                                  :: ncstatus

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Open external file

    ncstatus = nf90_open(path=trim(variable_filename),mode=                &
         & nf90_nowrite,ncid=ncfileid)

    ! Define local variables

    ncstatus     = nf90_inq_varid(ncfileid,'seaice',ncvarid)
    ncstatus     = nf90_get_var(ncfileid,ncvarid,workgrid_seaice)
    ncstatus     = nf90_inq_varid(ncfileid,'landmask',ncvarid)
    ncstatus     = nf90_get_var(ncfileid,ncvarid,workgrid_xland)

    ! Define local variable

    srcgrid%mask = real(workgrid_xland - workgrid_seaice)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !---------------------------------------------------------------------

    ! Loop through total number of destination grid coordinates

    do i = 1, dstgrid%ncoords

       ! Define local variable

       call kdtree2_n_nearest(tp=srcgrid%kdtree_grid,                      &
            & qv=dstgrid%grdloc(:,i),nn=1,results=sresults)

       ! Define local variable

       dstgrid%mask(i) = srcgrid%mask(sresults(1)%idx)

    end do ! do i = 1, dstgrid%ncoords

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_define_mask

  !=======================================================================

  ! interpolation_define_kdtree.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_kdtree(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Loop through total number of grid coordinates

    do j = 1, grid%ncoords

       ! Compute local variables

       grid%grdloc(1,j) = rearth_equator*cos(grid%xlat(j))*                &
            & cos(grid%xlong(j))
       grid%grdloc(2,j) = rearth_equator*cos(grid%xlat(j))*                &
            & sin(grid%xlong(j))
       grid%grdloc(3,j) = rearth_equator*sin(grid%xlat(j))

    end do ! do j = 1, grid%ncoords

    ! Initialize local variable

    grid%kdtree_grid => kdtree2_create(grid%grdloc,sort=.true.,            &
         & rearrange=.true.)

    !=====================================================================

  end subroutine interpolation_define_kdtree

  !=======================================================================

  ! interpolation_define_weights_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_weights_mpi(srcgrid,dstgrid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: srcgrid
    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    type(kdtree2_result),       dimension(srcgrid%neighbors)             :: sresults
    real(r_kind),               dimension(:,:,:),            allocatable :: mpi_weights
    real(r_kind),               dimension(:,:),              allocatable :: mpi_sresults_dis
    real(r_kind),               dimension(:,:),              allocatable :: mpi_scutoff
    real(r_kind)                                                         :: mpi_distthresh

    ! Define counting variables

    integer                                                              :: i, j, k, l 

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(mpi_weights))                                       &     
         & allocate(mpi_weights(dstgrid%ncoords,dstgrid%neighbors,         &
         & dstgrid%npasses))
    if(.not. allocated(mpi_scutoff))                                       &
         & allocate(mpi_scutoff(dstgrid%ncoords,dstgrid%npasses))
    if(.not. allocated(mpi_sresults_dis))                                  &
         & allocate(mpi_sresults_dis(dstgrid%ncoords,dstgrid%neighbors))

    !---------------------------------------------------------------------

    ! Initialize local variables

    dstgrid%weights = 0.0
    mpi_weights     = 0.0

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then
          
       ! Loop through local variable and proceed accordingly

       do i = 1, dstgrid%ncoords

          ! Define local variable
             
          call kdtree2_n_nearest(tp=srcgrid%kdtree_grid,                   &
               & qv=dstgrid%grdloc(:,i),nn=dstgrid%neighbors,              &
               & results=sresults)
             
          ! Define local variables
             
          mpi_distthresh                          =                        &
               & dstgrid%distthresh
          mpi_scutoff(i,1:dstgrid%npasses)        =                        &
               & dstgrid%scutoff(1:dstgrid%npasses)
          mpi_sresults_dis(i,1:dstgrid%neighbors) =                        &
               & sresults(1:dstgrid%neighbors)%dis
          dstgrid%grdnbors(i,1:dstgrid%neighbors) =                        &
               & sresults(1:dstgrid%neighbors)%idx

       end do ! do i = 1, dstgrid%ncoords
       
    end if ! if(mpi_procid .eq. mpi_masternode)
       
    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(mpi_distthresh,1,mpi_real,mpi_masternode,               &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(mpi_scutoff,(dstgrid%ncoords*dstgrid%npasses),          &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(mpi_sresults_dis,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%grdnbors,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_integer,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Check local variable and proceed accordingly

       if(dstgrid%mpi_count_begin(mpi_procid) .ne. 0 .and.                 &
            & dstgrid%mpi_count_end(mpi_procid) .ne. 0) then

          ! Loop through local variable and proceed accordingly
       
          do k = dstgrid%mpi_count_begin(mpi_procid),                      &
               & dstgrid%mpi_count_end(mpi_procid)
          
             ! Loop through total number of neighboring points on
             ! destination grid

             do i = 1, dstgrid%neighbors

                ! Loop through the total number of analysis passes to
                ! perform

                do j = 1, dstgrid%npasses

                   ! Compute local variable

                   mpi_weights(k,i,j) = exp((-1.0)*(                       &
                        & sqrt(mpi_sresults_dis(k,i))*                     &
                        & sqrt(mpi_sresults_dis(k,i)))/(4.0*               &
                        & mpi_scutoff(k,j)*mpi_distthresh*                 &
                        & mpi_distthresh))

                end do ! do j = 1, dstgrid%npasses

             end do ! do i = 1, dstgrid%neighbors

          end do ! do k = dstgrid%mpi_count_begin(mpi_procid),             &
                 ! dstgrid%mpi_count_end(mpi_procid)

       endif ! if(dstgrid%mpi_count_begin(mpi_procid) .ne. 0 .and.         &
             ! dstgrid%mpi_count_end(mpi_procid) .ne. 0)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_weights(1:dstgrid%ncoords,1:dstgrid%neighbors,1:   &
         & dstgrid%npasses),dstgrid%weights(1:dstgrid%ncoords,1:           &
         & dstgrid%neighbors,1:dstgrid%npasses),(dstgrid%ncoords*          &
         & dstgrid%neighbors*dstgrid%npasses),mpi_real,mpi_sum,            &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable on all compute tasks

    if(allocated(mpi_weights))      deallocate(mpi_weights)
    if(allocated(mpi_sresults_dis)) deallocate(mpi_sresults_dis)
    if(allocated(mpi_scutoff))      deallocate(mpi_scutoff)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_define_weights_mpi

  !=======================================================================

  ! interpolation_barnes_analysis_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_barnes_analysis_mpi(srcgrid,dstgrid,grid)

    ! Define variable passed to routine

    type(interpgrid)                                                     :: srcgrid
    type(grid_interface)                                                 :: grid

    ! Define variable returned by routine

    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    real(r_kind),               dimension(:),                allocatable :: mpi_dstgrid_var
    real(r_kind),               dimension(:),                allocatable :: workgrid
    real(r_kind)                                                         :: weights_sum

    ! Define counting variables

    integer                                                              :: i, j, k, l

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(mpi_dstgrid_var))                                   &
         & allocate(mpi_dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    mpi_dstgrid_var = 0.0

    ! Check local variable and proceed accordingly

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do k = dstgrid%mpi_count_begin(mpi_procid),                         &
            & dstgrid%mpi_count_end(mpi_procid)

          ! Loop through the total number of analysis passes
          
          do j = 1, dstgrid%npasses
             
             ! Initialize local variable
             
             weights_sum = 0.0
             
             ! Loop through total number of neighboring points on
             ! source grid

             do i = 1, dstgrid%neighbors
             
                ! Define local variable
                
                weights_sum = weights_sum + dstgrid%weights(k,i,j)
             
                ! Define local variable
             
                mpi_dstgrid_var(k) = mpi_dstgrid_var(k) +                  &
                     & (srcgrid_var(dstgrid%grdnbors(k,i)) -               &
                     & workgrid(k))*dstgrid%weights(k,i,j)
                
             end do ! do i = 1, dstgrid%neighbors

             ! Define local variable accordingly

             if(weights_sum .gt. barnes_weights_threshold .and. j .eq.     &
                  & 1) then

                ! Define local variable

                mpi_dstgrid_var(k) = mpi_dstgrid_var(k)/weights_sum

             else  ! if(weights_sum .gt. barnes_weights_threshold          &
                   ! .and. k .eq. 1)
                
                ! Define local variable

                mpi_dstgrid_var(k) = srcgrid_var(dstgrid%grdnbors(k,1))

             end if ! if(weights_sum .gt. barnes_weights_threshold         &
                    ! .and. k .eq. 1)
             
             ! Define local variable
             
             workgrid(k) = mpi_dstgrid_var(k)
                
          end do ! do j = 1, dstgrid%npasses

       end do ! do k = dstgrid%mpi_count_begin(mpi_procid),                &
              ! dstgrid%mpi_count_end(mpi_procid)

    endif ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid             &
          ! .le. dstgrid%mpi_maxprocid)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_dstgrid_var(1:dstgrid%ncoords),                    &
         & dstgrid_var(1:dstgrid%ncoords),dstgrid%ncoords,mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_dstgrid_var)) deallocate(mpi_dstgrid_var)
    if(allocated(workgrid))        deallocate(workgrid)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_barnes_analysis_mpi

  !=======================================================================

  ! interpolation_nearest_neighbor_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_nearest_neighbor_mpi(srcgrid,dstgrid)

    ! Define variable passed to routine

    type(interpgrid)                                                     :: srcgrid

    ! Define variable returned by routine

    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    type(kdtree2_result),       dimension(srcgrid%neighbors)             :: sresults
    real(r_kind),               dimension(:),                allocatable :: mpi_dstgrid_var
    integer,                    dimension(:,:),              allocatable :: mpi_srcgrid_idx

    ! Define counting variables

    integer                                                              :: i, j, k, l

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(mpi_srcgrid_idx))                                   &
         & allocate(mpi_srcgrid_idx(dstgrid%ncoords,dstgrid%neighbors))
    if(.not. allocated(mpi_dstgrid_var))                                   &
         & allocate(mpi_dstgrid_var(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    mpi_dstgrid_var = 0.0

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Loop through local variable and proceed accordingly

       do i = 1, dstgrid%ncoords

          ! Define local variable
             
          call kdtree2_n_nearest(tp=srcgrid%kdtree_grid,                   &
               & qv=dstgrid%grdloc(:,i),nn=dstgrid%neighbors,              &
               & results=sresults)
             
             ! Define local variable

          mpi_srcgrid_idx(i,1:dstgrid%neighbors) =                         &
               & sresults(1:dstgrid%neighbors)%idx

       end do ! i = 1, dstgrid%ncoords

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(mpi_srcgrid_idx,(dstgrid%ncoords*dstgrid%neighbors),    &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly
       
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do k = dstgrid%mpi_count_begin(mpi_procid),                         &
            & dstgrid%mpi_count_end(mpi_procid)

          ! Loop through total number of neighboring points on source
          ! grid

          do i = 1, dstgrid%neighbors

             ! Check local variable and proceed accordingly

             if(dstgrid%mask(k) .eq.                                       &
                  & srcgrid%mask(mpi_srcgrid_idx(k,i))) then

                ! Define local variable

                mpi_dstgrid_var(k) = srcgrid_var(mpi_srcgrid_idx(k,i))

                ! Exit loop

                goto 1000

             end if ! if(dstgrid%mask(k) .eq.                              &
                    ! srcgrid%mask(mpi_srcgrid_idx(k,i)))

          end do ! do i = 1, dstgrid%neighbors

          ! Define exit from loop

1000      continue

       end do ! do k = dstgrid%mpi_count_begin(mpi_procid),                &
              ! dstgrid%mpi_count_end(mpi_procid)

    endif ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.        &
          ! dstgrid%mpi_maxprocid)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_dstgrid_var(1:dstgrid%ncoords),                    &
         & dstgrid_var(1:dstgrid%ncoords),dstgrid%ncoords,mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_srcgrid_idx)) deallocate(mpi_srcgrid_idx)
    if(allocated(mpi_dstgrid_var)) deallocate(mpi_dstgrid_var)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_nearest_neighbor_mpi

  !=======================================================================

  ! interpolation_mpas_pressure_profile.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_mpas_pressure_profile(srcgrid,dstgrid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: srcgrid
    type(interpgrid)                                         :: dstgrid

    ! Define variables computed within routine

    real(r_double),             dimension(:,:),  allocatable :: workgrid_2d
    real(r_double),             dimension(:),    allocatable :: workgrid_1d
    real(r_kind),               dimension(:,:),  allocatable :: workgrid_pres
    real(r_kind),               dimension(:,:),  allocatable :: workgrid_pres_b
    real(r_kind),               dimension(:,:),  allocatable :: workgrid_pres_p
    real(r_kind),               dimension(:),    allocatable :: workgrid_psfc
    integer                                                  :: nvertcoords
    integer                                                  :: ncfileid
    integer                                                  :: ncvarid
    integer                                                  :: ncdimid
    integer                                                  :: ncstatus

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Open external file

    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define local variables

    ncstatus = nf90_inq_dimid(ncfileid,'nVertLevels',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,len=nvertcoords)

    ! Allocate memory for local variables

    if(.not. allocated(workgrid_2d))                                       &
         & allocate(workgrid_2d(nvertcoords,srcgrid%ncoords))
    if(.not. allocated(workgrid_1d))                                       &
         & allocate(workgrid_1d(srcgrid%ncoords))
    if(.not. allocated(workgrid_pres_b))                                   &
         & allocate(workgrid_pres_b(nvertcoords,srcgrid%ncoords))
    if(.not. allocated(workgrid_pres_p))                                   &
         & allocate(workgrid_pres_p(nvertcoords,srcgrid%ncoords))
    if(.not. allocated(workgrid_pres))                                     &
         & allocate(workgrid_pres(nvertcoords,srcgrid%ncoords))
    if(.not. allocated(workgrid_psfc))                                     &
         & allocate(workgrid_psfc(srcgrid%ncoords))

    ! Define local variables

    ncstatus        = nf90_inq_varid(ncfileid,'pressure_base',ncvarid)
    ncstatus        = nf90_get_var(ncfileid,ncvarid,workgrid_2d)
    workgrid_pres_b = real(workgrid_2d)
    ncstatus        = nf90_inq_varid(ncfileid,'pressure_p',ncvarid)
    ncstatus        = nf90_get_var(ncfileid,ncvarid,workgrid_2d)
    workgrid_pres_p = real(workgrid_2d)

    ! Compute local variable

    workgrid_pres = workgrid_pres_b + workgrid_pres_p

    ! Define local variable

    srcgrid%pressure_profile = workgrid_pres

    ! Define local variables

    ncstatus                 =                                             &
         & nf90_inq_varid(ncfileid,'surface_pressure',ncvarid)
    ncstatus                 =                                             &
         & nf90_get_var(ncfileid,ncvarid,workgrid_1d)
    workgrid_psfc            = real(workgrid_1d)
    srcgrid%surface_pressure = workgrid_psfc

    ! Deallocate memory for local variables

    if(allocated(workgrid_2d))     deallocate(workgrid_2d)
    if(allocated(workgrid_1d))     deallocate(workgrid_1d)
    if(allocated(workgrid_pres_b)) deallocate(workgrid_pres_b)
    if(allocated(workgrid_pres_p)) deallocate(workgrid_pres_p)
    if(allocated(workgrid_pres))   deallocate(workgrid_pres)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_mpas_pressure_profile

  !=======================================================================

  ! interpolation_mpas_height_profile.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_mpas_height_profile(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    ! Define variables computed within routine

    real(r_double),             dimension(:,:),  allocatable :: workgrid
    integer                                                  :: nvertcoords
    integer                                                  :: ncfileid
    integer                                                  :: ncvarid
    integer                                                  :: ncdimid
    integer                                                  :: ncstatus

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Open external file

    ncstatus = nf90_open(path=trim(variable_filename),                     &
         & mode=nf90_nowrite,ncid=ncfileid)

    ! Define local variables

    ncstatus = nf90_inq_dimid(ncfileid,'nVertLevelsP1',ncdimid)
    ncstatus = nf90_inquire_dimension(ncfileid,ncdimid,len=nvertcoords)

    ! Allocate memory for local variable

    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(nvertcoords,grid%ncoords))

    ! Define local variables

    ncstatus = nf90_inq_varid(ncfileid,'zgrid',ncvarid)
    ncstatus = nf90_get_var(ncfileid,ncvarid,workgrid)
    ncstatus = nf90_close(ncfileid)

    ! Define local variable

    grid%surface_height(1:grid%ncoords) =                                  &
         & real(workgrid(1,1:grid%ncoords))

    ! Loop though vertical coordinate

    do k = 1, grid%nvertlevs

       ! Compute local variable

       grid%height_profile(k,:) = real(workgrid(k,:) + workgrid(k+1,:))    &
            & /2.0

    end do ! do k = 1, grid%nvertlevs

    ! Deallocate memory for local variables

    if(allocated(workgrid)) deallocate(workgrid)    

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_mpas_height_profile

  !=======================================================================

  ! interpolation_mpas_thetae_profile.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_mpas_thetae_profile(srcgrid)

    ! Define variables passed to routine

    type(interpgrid)                                                      :: srcgrid

    ! Define variables computed within routine

    real(r_double),          dimension(srcgrid%nvertlevs,srcgrid%ncoords) :: pressure
    real(r_double),          dimension(srcgrid%nvertlevs,srcgrid%ncoords) :: theta
    real(r_double),          dimension(srcgrid%nvertlevs,srcgrid%ncoords) :: qvapor
    real(r_double),          dimension(srcgrid%nvertlevs,srcgrid%ncoords) :: temperature
    real(r_double),          dimension(srcgrid%nvertlevs,srcgrid%ncoords) :: workgrid
    real(r_double)                                                        :: vappres
    real(r_double)                                                        :: tlcl
    integer                                                               :: ncfileid
    integer                                                               :: ncvarid
    integer                                                               :: ncdimid
    integer                                                               :: ncstatus

    ! Define counting variables

    integer                                                               :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

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

    ! Loop through all vertical coordinates and proceed accordingly

    do k = 1, srcgrid%nvertlevs

       ! Loop through all horizontal coordinates and proceed
       ! accordingly

       do i = 1, srcgrid%ncoords

          ! Compute local variables

          vappres = pressure(k,i)*qvapor(k,i)/(dble(eps)+qvapor(k,i))/     &
               & dble(100.0)
          tlcl    = 55.0 + 2840.0/(3.5*log(temperature(k,i)) -             &
               & log(vappres) - 4.805)
          
          ! Rescale local variable accordingly

          tlcl = min(tlcl,temperature(k,i))

          ! Compute local variable

          srcgrid%thetae_profile(k,i) = theta(k,i)*                        &
               & exp(((3376./tlcl)-2.54)*qvapor(k,i)*(1.0 +                &
               & 0.81*qvapor(k,i)))

       end do ! do i = 1, srcgrid%ncoords

    end do ! do k = 1, srcgrid%nvertlevs

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_mpas_thetae_profile

  !=======================================================================

  ! interpolation_interppres_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_interppres_mpi(interppres_grid,dstgrid)

    ! Define variables passed to routine

    type(interppres)                                         :: interppres_grid
    type(interpgrid)                                         :: dstgrid

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),  allocatable :: mpi_profile_var
    
    ! Define counting variables

    integer                                                  :: i, j, k, l

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%nvertlevs,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(interppres_grid%nsig,1,mpi_integer,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%pressure_profile,                               &
         & (dstgrid%ncoords*dstgrid%nvertlevs),mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(interppres_grid%profile_pres,                           &
         & (dstgrid%ncoords*interppres_grid%nsig),mpi_real,mpi_masternode, &
         & mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(mpi_profile_var))                                   &
         & allocate(mpi_profile_var(dstgrid%ncoords,interppres_grid%nsig))
    
    ! Initialize local variable

    mpi_profile_var = 0.0
    
    !---------------------------------------------------------------------

    ! Check local variable and proceed accordingly
    
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do i = dstgrid%mpi_count_begin(mpi_procid),                         &
            & dstgrid%mpi_count_end(mpi_procid)

          ! Compute local variable
          
          call calculate_profile_linearlogpinterpolation(                  &
               & dstgrid%nvertlevs,interppres_grid%nsig,                   &
               & dstgrid%pressure_profile(1:dstgrid%nvertlevs,i),          &
               & dstgrid%anlysvar(i,1:dstgrid%nvertlevs),                  &
               & interppres_grid%profile_pres(i,1:                         &
               & interppres_grid%nsig),mpi_profile_var(i,                  &
               & 1:interppres_grid%nsig))

          ! Loop through vertical coordinate
          
          do k = 1, interppres_grid%nsig
             
             ! Update local variable accordingly
             
             if(dstgrid%surface_pressure(i) .lt.                           &
                  & interppres_grid%profile_pres(i,k)) then

                ! Update local variable
                
                mpi_profile_var(i,k) = 1.e30
                   
             end if ! if(dstgrid%surface_pressure(i) .lt.                  &
                    ! dstgrid%pressure_profile(i,k))
             
          end do ! do k = 1, interppres_grid%nsig

       end do ! do i = dstgrid%mpi_count_begin(mpi_procid),                &
              ! dstgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.       &
           ! dstgrid%mpi_maxprocid)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_profile_var(1:dstgrid%ncoords,                     &
         & 1:interppres_grid%nsig),interppres_grid%profile_var(            &
         & 1:dstgrid%ncoords,1:interppres_grid%nsig),(dstgrid%ncoords*     &
         & interppres_grid%nsig),mpi_real,mpi_sum,mpi_masternode,          &
         & mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then

       ! Loop through vertical coordinate

       do k = 1, interppres_grid%nsig

          ! Print message to user

          if(debug) write(6,500) k,                                        &
               & minval(interppres_grid%profile_pres(:,k)),                &
               & maxval(interppres_grid%profile_pres(:,k)),                & 
               & minval(interppres_grid%profile_var(:,k)),                 &
               & maxval(interppres_grid%profile_var(:,k),                  &
               & interppres_grid%profile_var(:,k) .ne. 1.e30)

       end do ! do k = 1, interppres_grid%nsig
 
    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(mpi_profile_var))  deallocate(mpi_profile_var)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_INTERPPRES_MPI: (level, min/max pressure, ',     &
         & 'min/max variable)', i6, f, '/', f, f, '/', f)

    !=====================================================================

  end subroutine interpolation_interppres_mpi

  !=======================================================================

  ! interpolation_interphght_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_interphght_mpi(interphght_grid,grid)

    ! Define variables passed to routine

    type(interphght)                                         :: interphght_grid
    type(interpgrid)                                         :: grid

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),  allocatable :: mpi_profile_var
    
    ! Define counting variables

    integer                                                  :: i, j, k, l

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%ncoords,1,mpi_integer,mpi_masternode,              &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(interphght_grid%ncoords,1,mpi_integer,mpi_masternode,   &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%nvertlevs,1,mpi_integer,mpi_masternode,            &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(interphght_grid%nsig,1,mpi_integer,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%anlysvar,                                          &
         & (grid%ncoords*grid%nvertlevs),mpi_real,                         &         
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%height_profile,                                    &
         & (grid%ncoords*grid%nvertlevs),mpi_real,                         &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(interphght_grid%profile_hght,                           &
         & (interphght_grid%ncoords*interphght_grid%nsig),mpi_real,        &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(mpi_profile_var))                                   &
         & allocate(mpi_profile_var(grid%ncoords,interphght_grid%nsig))
 
    ! Initialize local variable

    mpi_profile_var = 0.0

    !---------------------------------------------------------------------

    ! Check local variable and proceed accordingly
    
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & grid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do i = grid%mpi_count_begin(mpi_procid),                            &
            & grid%mpi_count_end(mpi_procid)

          ! Define local variable

          interphght_grid%profile_hght(i,1:interphght_grid%nsig) =         &
               & interphght_grid%profile_hght(i,1:interphght_grid%nsig) +  &
               & grid%surface_height(i)

          ! Compute local variable

          call calculate_profile_linearinterpolation(                      &
               & grid%nvertlevs,interphght_grid%nsig,                      &
               & grid%height_profile(1:grid%nvertlevs,i),                  &
               & grid%anlysvar(i,1:grid%nvertlevs),                        &
               & interphght_grid%profile_hght(i,                           &
               & 1:interphght_grid%nsig),mpi_profile_var(i,                &
               & 1:interphght_grid%nsig))

       end do ! do i = grid%mpi_count_begin(mpi_procid),                   &
              ! grid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.       &
           ! grid%mpi_maxprocid)

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_profile_var(1:grid%ncoords,                        &
         & 1:interphght_grid%nsig),interphght_grid%profile_var(            &
         & 1:grid%ncoords,1:interphght_grid%nsig),                         &
         & (grid%ncoords*interphght_grid%nsig),mpi_real,mpi_sum,           &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then

       ! Loop through vertical coordinate

       do k = 1, interphght_grid%nsig

          ! Print message to user

          if(debug) write(6,503) k,                                        &
               & minval(interphght_grid%profile_hght(:,k)),                &
               & maxval(interphght_grid%profile_hght(:,k)),                & 
               & minval(interphght_grid%profile_var(:,k)),                 &
               & maxval(interphght_grid%profile_var(:,k),                  &
               & interphght_grid%profile_var(:,k) .ne. 1.e30)

       end do ! do k = 1, interphght_grid%nsig
 
    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(mpi_profile_var))  deallocate(mpi_profile_var)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_INTERPHGHT_MPI: A total of ', i6, ' loops will', &
         & ' be required to accomplish tasks.')
501 format('INTERPOLATION_INTERPHGHT_MPI: Sending grid location ', i9,     &
         & ' of ', i9, ' to compute task.')
502 format('INTERPOLATION_INTERPHGHT_MPI: Receiving grid location ', i9,   &
         & ' of ', i9, ' from compute task.')
503 format('INTERPOLATION_INTERPHGHT_MPI: (level, min/max height, ',       &
         & 'min/max variable)', i6, f, '/', f, f, '/', f)

    !=====================================================================

  end subroutine interpolation_interphght_mpi

  !=======================================================================

  ! interpolation_interpthetae_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_interpthetae_mpi(interpthee_grid,grid)

    ! Define variables passed to routine

    type(interpthetae)                                       :: interpthee_grid
    type(interpgrid)                                         :: grid

    ! Define variables computed within routine

    real(r_kind),               dimension(:,:),  allocatable :: mpi_profile_var
    
    ! Define counting variables

    integer                                                  :: i, j, k, l

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()
    call init_constants(.true.)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%ncoords,1,mpi_integer,mpi_masternode,              &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(interpthee_grid%ncoords,1,mpi_integer,mpi_masternode,   &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%nvertlevs,1,mpi_integer,mpi_masternode,            &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(interpthee_grid%nsig,1,mpi_integer,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%anlysvar,                                          &
         & (grid%ncoords*grid%nvertlevs),mpi_real,                         &         
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%height_profile,                                    &
         & (grid%ncoords*grid%nvertlevs),mpi_real,                         &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(interpthee_grid%profile_thetae,                         &
         & (interpthee_grid%ncoords*interpthee_grid%nsig),mpi_real,        &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(mpi_profile_var))                                   &
         & allocate(mpi_profile_var(grid%ncoords,interpthee_grid%nsig))
 
    ! Initialize local variable

    mpi_profile_var = 0.0

    !---------------------------------------------------------------------

    ! Check local variable and proceed accordingly
    
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & grid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do i = grid%mpi_count_begin(mpi_procid),                            &
            & grid%mpi_count_end(mpi_procid)


          ! Compute local variable

          call calculate_profile_linearinterpolation(                      &
               & grid%nvertlevs,interpthee_grid%nsig,                      &
               & grid%thetae_profile(1:grid%nvertlevs,i),                  &
               & grid%anlysvar(i,1:grid%nvertlevs),                        &
               & interpthee_grid%profile_thetae(i,                         &
               & 1:interpthee_grid%nsig),mpi_profile_var(i,                &
               & 1:interpthee_grid%nsig))

       end do ! do i = grid%mpi_count_begin(mpi_procid),                   &
              ! grid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.       &
           ! grid%mpi_maxprocid)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_profile_var(1:grid%ncoords,                        &
         & 1:interpthee_grid%nsig),interpthee_grid%profile_var(            &
         & 1:grid%ncoords,1:interpthee_grid%nsig),                         &
         & (grid%ncoords*interpthee_grid%nsig),mpi_real,mpi_sum,           &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then

       ! Loop through vertical coordinate

       do k = 1, interpthee_grid%nsig

          ! Print message to user

          if(debug) write(6,503) k,                                        &
               & minval(interpthee_grid%profile_thetae(:,k)),              &
               & maxval(interpthee_grid%profile_thetae(:,k)),              & 
               & minval(interpthee_grid%profile_var(:,k)),                 &
               & maxval(interpthee_grid%profile_var(:,k),                  &
               & interpthee_grid%profile_var(:,k) .ne. 1.e30)

       end do ! do k = 1, interpthee_grid%nsig
 
    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(mpi_profile_var))  deallocate(mpi_profile_var)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_INTERPTHETAE_MPI: A total of ', i6, ' loops ',   &
         & 'will be required to accomplish tasks.')
501 format('INTERPOLATION_INTERPTHETAE_MPI: Sending grid location ', i9,   &
         & ' of ', i9, ' to compute task.')
502 format('INTERPOLATION_INTERPTHETAE_MPI: Receiving grid location ', i9, &
         & ' of ', i9, ' from compute task.')
503 format('INTERPOLATION_INTERPTHETAE_MPI: (level, min/max height, ',     &
         & 'min/max variable)', i6, f, '/', f, f, '/', f)

    !=====================================================================

  end subroutine interpolation_interpthetae_mpi

  !=======================================================================

  ! interpolation_define_kdtree_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_kdtree_mpi(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),   allocatable :: mpi_grdloc

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%xlong,grid%ncoords,mpi_real,mpi_masternode,        &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%xlat,grid%ncoords,mpi_real,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Allocate memory for local variable

    if(.not. allocated(mpi_grdloc)) allocate(mpi_grdloc(3,grid%ncoords))

    ! Initialize local variables

    mpi_grdloc  = 0.0
    grid%grdloc = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Check local variable and proceed accordingly

       if(grid%mpi_count_begin(mpi_procid) .ne. 0 .and.                    &
            & grid%mpi_count_end(mpi_procid) .ne. 0) then

          ! Loop through local variable and proceed accordingly
       
          do k = grid%mpi_count_begin(mpi_procid),                         &
               & grid%mpi_count_end(mpi_procid)
          
             ! Compute local variables

             mpi_grdloc(1,k) = rearth_equator*cos(grid%xlat(k))*           &
                  & cos(grid%xlong(k))
             mpi_grdloc(2,k) = rearth_equator*cos(grid%xlat(k))*           &
                  & sin(grid%xlong(k))
             mpi_grdloc(3,k) = rearth_equator*sin(grid%xlat(k))

          end do ! do k = grid%mpi_count_begin(mpi_procid),                &
                 ! grid%mpi_count_end(mpi_procid)

       end if ! if(grid%mpi_count_begin(mpi_procid) .ne. 0 .and.           &
              ! grid%mpi_count_end(mpi_procid) .ne. 0)

    end if ! if(mpi_procid .ne. mpi_masternode)

    ! Define local variable

    call mpi_reduce(mpi_grdloc(1:3,1:grid%ncoords),grid%grdloc(1:3,1:      &
         & grid%ncoords),(3*grid%ncoords),mpi_real,mpi_sum,                &
         & mpi_masternode,mpi_comm_world,mpi_ierror)  

    ! Deallocate memory for local variable

    if(allocated(mpi_grdloc)) deallocate(mpi_grdloc)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize local variable

       grid%kdtree_grid => kdtree2_create(grid%grdloc,sort=.true.,         &
            & rearrange=.true.)

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine interpolation_define_kdtree_mpi

  !=======================================================================

  ! interpolation_smoothanalysis.f90:

  !-----------------------------------------------------------------------
  
  subroutine interpolation_smoothanalysis(array)

    ! Define variable passed to routine

    real(r_kind),               dimension(nlon*nlat)         :: array

    ! Define variables computed within routine

    real(r_kind),               dimension(nlon,nlat)         :: workgrid
    real(r_kind),               dimension(nlon*nlat)         :: smooth

    ! Define counting variables

    integer                                                  :: i, j, k
    integer                                                  :: count

    !=====================================================================

    ! Compute local variable

    call calculate_spectralsmoothingfunction(smooth)

    !---------------------------------------------------------------------

    ! Initialize local variable

    count = 1
    
    ! Loop through meridional horizontal coordinate
    
    do j = 1, nlat
       
       ! Loop through zonal horizontal coordinate
       
       do i = 1, nlon

          ! Define local variable

          workgrid(i,j) = real(array(count))

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, nlon

    end do ! do j = 1, nlat

    ! Compute local variable

    call specsmooth(workgrid,smooth,'GAU')

    ! Initialize local variable

    count = 1
    
    ! Loop through meridional horizontal coordinate
    
    do j = 1, nlat
       
       ! Loop through zonal horizontal coordinate
       
       do i = 1, nlon

          ! Define local variable

          array(count) = workgrid(i,j)

          ! Update counting variable

          count = count + 1

       end do ! do i = 1, nlon

    end do ! do j = 1, nlat

    !=====================================================================

  end subroutine interpolation_smoothanalysis

  !=======================================================================

  ! define_scaling_coefficients.f90:

  !-----------------------------------------------------------------------

  subroutine define_scaling_coefficients(grid)

    ! Define variable passed to routine

    type(interpgrid)                                         :: grid

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variable

    grid%scutoff(1) = 1.0

    ! Loop through total number of threshold cutoff values

    do k = 1, grid%npasses

       ! Compute local variable

       grid%scutoff(k) = grid%scutoff(1)/(10**(real(k-1)))

    end do ! do k = 1, grid%npasses

    !=====================================================================

  end subroutine define_scaling_coefficients

  !=======================================================================

  ! calculate_profile_linearlogpinterpolation.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_profile_linearlogpinterpolation(src_nlev,dst_nlev,  &
       & src_pres,src_var,dst_pres,dst_var)

    ! Define array dimension variables

    integer                                                  :: src_nlev
    integer                                                  :: dst_nlev

    ! Define variables passed to routine

    real(r_kind),              dimension(src_nlev)           :: src_pres
    real(r_kind),              dimension(src_nlev)           :: src_var
    real(r_kind),              dimension(dst_nlev)           :: dst_pres

    ! Define variables returned by routine

    real(r_kind),              dimension(dst_nlev)           :: dst_var

    ! Define variables computed within routine

    integer                                                  :: index_bottom
    integer                                                  :: index_top

    ! Define counting variables

    integer                                                  :: i, j, k, l

    !=====================================================================

    ! Loop through vertical coordinate

    do l = 1, dst_nlev

       ! Initialize local variables

       index_bottom = 0
       index_top    = 0

       ! Loop through vertical coordinate
    
       do k = 2, src_nlev 

          ! Define local variables accordingly

          if(dst_pres(l) .le. src_pres(k-1) .and. dst_pres(l) .gt.         &
               & src_pres(k)) then

             ! Define local variables

             index_bottom = k - 1
             index_top    = k

             ! Exit loop

             goto 100

          end if ! if(dst_pres(l) .le. src_pres(k-1) .and. dst_pres(l)     &
                 ! .gt. src_pres(k))

       end do ! do k = 2, src_nlev 

       ! Define exit from loop

100    continue

       ! Update local variables accordingly

       if(index_bottom .eq. 0 .or. index_top .eq. 0) then

          ! Update local variables

          index_bottom = 1
          index_top    = src_nlev

       end if ! if(index_bottom .eq. 0 .or. index_top .eq. 0)

       ! Compute local variable

       dst_var(l) = src_var(index_bottom) + (src_var(index_top) -          & 
            & src_var(index_bottom))*((log(dst_pres(l)) -                  &
            & log(src_pres(index_bottom)))/(log(src_pres(index_top)) -     &
            & log(src_pres(index_bottom))))

    end do ! do k = 1, src_nlev

    !=====================================================================

  end subroutine calculate_profile_linearlogpinterpolation

  !=======================================================================

  ! calculate_profile_linearinterpolation.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_profile_linearinterpolation(src_nlev,dst_nlev,      &
       & src_hght,src_var,dst_hght,dst_var)

    ! Define array dimension variables

    integer                                                  :: src_nlev
    integer                                                  :: dst_nlev

    ! Define variables passed to routine

    real(r_kind),              dimension(src_nlev)           :: src_hght
    real(r_kind),              dimension(src_nlev)           :: src_var
    real(r_kind),              dimension(dst_nlev)           :: dst_hght

    ! Define variables returned by routine

    real(r_kind),              dimension(dst_nlev)           :: dst_var

    ! Define variables computed within routine

    logical                                                  :: interp
    real(r_kind)                                             :: w1
    real(r_kind)                                             :: w2
    integer                                                  :: src_idx

    ! Define counting variables

    integer                                                  :: i, j, k, l

    !=====================================================================

    ! Loop through vertical coordinate

    do k = 1, dst_nlev

       ! Initialize local variables

       src_idx = 1
       interp  = .false.

       ! Loop through local variable and proceed accordingly

       do while((.not. interp) .and. src_idx .lt. src_nlev)

          ! Define local variables accordingly

          if(src_hght(src_idx) .le. dst_hght(k) .and. src_hght(src_idx+1)  &
               & .ge. dst_hght(k)) then

             ! Define local variable
             
             interp = .true. 

             ! Compute local variables

             w1         = (src_hght(src_idx+1) - dst_hght(k))/             &
                  & (src_hght(src_idx+1) - src_hght(src_idx))
             w2         = 1.0 - w1
             dst_var(k) = w1*src_var(src_idx) + w2*src_var(src_idx+1)

          end if ! if(src_hght(src_idx) .le. dst_hght(k)                   &
                 ! .and. src(src_idx + 1) .ge. dst_hght(k))
             
          ! Update local variable

          src_idx = src_idx + 1

       end do ! do while((.not. interp) .and. src_idx .lt. src_nlev)

       ! Check local variable and proceed accordingly

       if(.not. interp) then

          ! Compute local variables

          w1         = (src_hght(src_nlev) - dst_hght(k))/                 &
               & (src_hght(src_nlev) - src_hght(1))
          w2         = 1.0 - w1
          dst_var(k) = w1*src_var(1) + w2*src_var(src_nlev)

       end if ! if(.not. interp)

    end do ! do k = 1, dst_nlev

    !=====================================================================

  end subroutine calculate_profile_linearinterpolation

  !=======================================================================

  ! calculate_spectralsmoothingfunction.f90:  

  !-----------------------------------------------------------------------

  subroutine calculate_spectralsmoothingfunction(smooth)

    ! Define variables returned by routine

    real(r_kind),              dimension(nlon*nlat)                   :: smooth

    ! Define counting variables

    integer                                                           :: i, j, k

    !=====================================================================

    ! Initialize local variable

    smooth = 0.0

    ! Loop through total triangular trunction and proceed accordingly

    do k = 1, (grid_mtrunc + 1)

       ! Compute local variable

       smooth(k) = exp(-(float(k-1)/real(grid_mtrunc-k+1)**2.0))

    end do ! do k = 1, (grid_mtrunc + 1)

    !=====================================================================

    ! Return calculated variables

    return

    !=====================================================================

  end subroutine calculate_spectralsmoothingfunction

  !=======================================================================

end module interpolation_interface
