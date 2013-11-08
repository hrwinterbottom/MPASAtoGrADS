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

module grads_interface

  !=======================================================================
  
  ! Define associated modules and subroutines
  
  !-----------------------------------------------------------------------
  
  use constants
  use kinds

  !-----------------------------------------------------------------------

  use namelist
  use netcdf
  use variable_interface

  !-----------------------------------------------------------------------
  
  implicit none

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! grads_interface_define_timestamp.f90:

  !-----------------------------------------------------------------------

  subroutine grads_interface_define_timestamp(grid_time,timestamp)

    ! Define variables passed to subroutine

    character(len=19)                         :: grid_time

    ! Define variables returned by subroutine

    character(len=15)                         :: timestamp

    ! Define variables computed within subroutine

    character(len=3)                          :: timestamp_month
    character(len=4)                          :: year
    character(len=2)                          :: month
    character(len=2)                          :: day
    character(len=2)                          :: hour
    character(len=2)                          :: minute

    !=================================================================

    ! Define local variables

    year   = grid_time(1:4)
    month  = grid_time(6:7)
    day    = grid_time(9:10)
    hour   = grid_time(12:13)
    minute = grid_time(15:16)

    ! Define local variable accordingly

    if(month .eq. '01') timestamp_month = 'JAN'
    if(month .eq. '02') timestamp_month = 'FEB'
    if(month .eq. '03') timestamp_month = 'MAR'
    if(month .eq. '04') timestamp_month = 'APR'
    if(month .eq. '05') timestamp_month = 'MAY'
    if(month .eq. '06') timestamp_month = 'JUN'
    if(month .eq. '07') timestamp_month = 'JUL'
    if(month .eq. '08') timestamp_month = 'AUG'
    if(month .eq. '09') timestamp_month = 'SEP'
    if(month .eq. '10') timestamp_month = 'OCT'
    if(month .eq. '11') timestamp_month = 'NOV'
    if(month .eq. '12') timestamp_month = 'DEC'

    ! Compute local variable

    write(timestamp,500) hour, minute, day, timestamp_month, year

    !=================================================================

    ! Return calculated values

    return

    !=================================================================

    ! Define format statements

500 format(a2,":",a2,"Z",a2,a3,a4)

    !=================================================================

  end subroutine grads_interface_define_timestamp

  !=======================================================================

  ! grads_interface_descriptor.f90:

  !-----------------------------------------------------------------------

  subroutine grads_interface_descriptor(grid)

    ! Define variables passed to routine

    type(variable_info),     dimension(number_of_var)        :: grid

    ! Define variables computed within routine

    character(len=50)                                        :: fmt
    character(len=19)                                        :: grid_time
    character(len=15)                                        :: grads_timestamp
    character(len=8)                                         :: grid_output_interval
    real(r_double),            dimension(:,:),   allocatable :: grid_variable_real
    real(r_kind),              dimension(:),     allocatable :: dstgrid_xlat
    real(r_kind),              dimension(:),     allocatable :: dstgrid_slat
    real(r_kind),              dimension(:),     allocatable :: dstgrid_wlat
    real(r_kind)                                             :: lon_min
    integer                                                  :: grid_variable_integer
    integer                                                  :: ncfileid
    integer                                                  :: ncvarid
    integer                                                  :: ncdimid
    integer                                                  :: ncstatus

    ! Define counting variables

    integer                                                  :: i, j, k

    !===================================================================== 

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Open external file

    ncstatus = nf90_open(path=trim(variable_filename),mode=nf90_nowrite,   &
         & ncid=ncfileid)

    ! Define local variable

    ncstatus = nf90_inq_varid(ncfileid,'xtime',ncvarid)

    ! Ingest local variable

    ncstatus = nf90_get_var(ncfileid,ncvarid,grid_time)

    ! Close external file

    ncstatus = nf90_close(ncfileid)

    ! Define local variable

    call grads_interface_define_timestamp(grid_time,grads_timestamp)

    ! Open external file

    open(999,file='mpastograds.ctl',form='formatted')

    ! Write values to external file

    write(999,500) 'mpastograds.bin'
    write(999,501)
    write(999,502)

    ! Allocate memory for local variables

    if(.not. allocated(dstgrid_xlat)) allocate(dstgrid_xlat(nlat))
    if(.not. allocated(dstgrid_slat)) allocate(dstgrid_slat(nlat))
    if(.not. allocated(dstgrid_wlat)) allocate(dstgrid_wlat(nlat))

    ! Define local variables

    call gausslat(nlat,dstgrid_slat,dstgrid_wlat)

    ! Compute local variablde

    dstgrid_xlat = (acos(dstgrid_slat) - pi/2.0)*rad2deg

    ! Write values to external file

    write(999,503) nlon, rlon_min, dx
    write(999,504) nlat
    write(999,505) dstgrid_xlat

    ! Deallocate memory for local variables

    if(allocated(dstgrid_xlat)) deallocate(dstgrid_xlat)
    if(allocated(dstgrid_slat)) deallocate(dstgrid_slat)
    if(allocated(dstgrid_wlat)) deallocate(dstgrid_wlat)

    ! Define local variable

    write(fmt,'("("i"(f13.5))")') nlev

    ! Write values to external file

    if(.not. is_pinterp .and. .not. is_zinterp .and. .not. is_teinterp)   &   
         & write(999,506) nlev
    if(is_pinterp .or. is_zinterp .or. is_teinterp .and. nlev .gt. 1)     &
         & write(999,507) nlev

    ! Check local variable and proceed accordingly

    if(nlev .gt. 1) then

       ! Write values to external file

       if(is_pinterp)  write(999,(adjustl(fmt)))                          &
            & (plevs(k)/100.0,k = 1, nlev)
       if(is_zinterp)  write(999,(adjustl(fmt)))                          &
            & (zlevs(k),k = 1, nlev)
       if(is_teinterp) write(999,(adjustl(fmt)))                          &
            & (televs(k),k = 1, nlev)

    else   ! if(nlev .gt. 1) 

       ! Write values to external file

       if(is_pinterp)  write(999,513) plevs(1)/100.0
       if(is_zinterp)  write(999,513) zlevs(1)
       if(is_teinterp) write(999,513) televs(1)

    end if ! if(nlev .gt. 1) 

    ! Write values to external file

    write(999,509) grads_timestamp

    ! Write values to external file

    write(999,510) number_of_var
    do i = 1, number_of_var
       write(999,511) grid(i)%grads_id, grid(i)%zdim, 99,                  &
            & grid(i)%grads_string
    end do ! do i = 1, number_of_var
    write(999,512)

    ! Close external file

    close(999)

    !---------------------------------------------------------------------

    ! Define format statements

500 format('dset ', a)
501 format('undef 1.e30')
502 format('options yrev')
503 format('xdef ', i6, ' linear ', f13.5,1x,f13.5)
504 format('ydef ', i6, ' levels ')
505 format(f13.5)
506 format('zdef ', i6, ' linear   0.0  1.0')
507 format('zdef ', i6, ' levels')
508 format(f13.5)
509 format('tdef 1 linear ', a15, ' 1hr')
510 format('vars ', i4)
511 format(a,1x,i2,1x,i6,1x,a)
512 format('endvars')
513 format('zdef 1 levels ',f13.5)

    !===================================================================== 

  end subroutine grads_interface_descriptor

  !=======================================================================

  ! grads_interface_binary.f90:

  !-----------------------------------------------------------------------

  subroutine grads_interface_binary(grid)

    ! Define variables passed to routine

    type(variable_info),     dimension(number_of_var)        :: grid

    ! Define counting variables

    integer                                                  :: i, j, k

    !===================================================================== 
 
    ! Open external file

    open(999,file='mpastograds.bin',form='unformatted',status='unknown',   &
         & recordtype='stream')

    ! Loop through each variable and proceed accordingly

    do i = 1, number_of_var

       ! Write variable to external file

       write(999) grid(i)%array
    
       ! Deallocate memory for local variable

       if(allocated(grid(i)%array)) deallocate(grid(i)%array)

    end do !  do i = 1, number_of_var

    ! Close external file

    close(999)

    !===================================================================== 

  end subroutine grads_interface_binary

  !=======================================================================

end module grads_interface
